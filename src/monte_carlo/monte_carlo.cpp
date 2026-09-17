#include "monte_carlo/monte_carlo.hpp"
#include <fstream>
#include <random>
#include <limits>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
namespace lv::monte_carlo {
std::vector<RunResult> run(const dynamics::VehicleParameters& base,const sensors::SensorConfig& sc,const Config& c,
                           std::vector<ConsistencySample>* consistency_samples) {
    std::vector<RunResult> out(c.runs);
    std::vector<std::vector<ConsistencySample>> traces(
        consistency_samples ? static_cast<std::size_t>(c.runs) : 0U);
#pragma omp parallel for if(c.threads > 1)
    for(int i=0;i<c.runs;++i) {
        std::mt19937_64 rng(c.seed+static_cast<std::uint64_t>(i)); std::normal_distribution<double> w(0,1);
        double wind_scale=1.0+c.wind_sigma*w(rng), thrust_scale=1.0+c.thrust_sigma*w(rng), dry_mass_scale=1.0+c.dry_mass_sigma*w(rng), aero_scale=1.0+c.aero_sigma*w(rng), bias_scale=1.0+c.sensor_bias_sigma*w(rng);
        double misalignment_deg=c.thrust_misalignment_sigma_deg*w(rng), cg_offset_m=c.cg_offset_sigma_m*w(rng);
        if(c.navigation_only) {
            // Deterministic IMU propagation diagnostic: remove plant
            // dispersions and sensor noise from this special mode.
            wind_scale=1.0; thrust_scale=1.0; dry_mass_scale=1.0;
            aero_scale=1.0; bias_scale=0.0; misalignment_deg=0.0; cg_offset_m=0.0;
        }
        auto p=base; p.dry_mass_kg*=std::max(0.1,dry_mass_scale);
        constexpr double deg2rad=0.017453292519943295; p.thrust_misalignment_pitch_rad=misalignment_deg*deg2rad; p.thrust_misalignment_yaw_rad=misalignment_deg*deg2rad;
        // A transverse CG offset creates a physical thrust moment; an axial
        // offset would be collinear with nominal thrust and create no moment.
        p.thrust_application_point_b_m.y()+=cg_offset_m;
        for(auto& q:p.aero.drag_vs_mach) q.value*=aero_scale;
        for(auto& q:p.wind) q.velocity_i_mps*=wind_scale;
        for(auto& m:p.motor) m.vacuum_thrust_n*=thrust_scale;
        auto sensor_config=sc;
        if(c.navigation_only) {
            sensor_config.accel_bias_sigma.setZero(); sensor_config.gyro_bias_sigma.setZero();
            sensor_config.accel_noise_sigma.setZero(); sensor_config.gyro_noise_sigma.setZero();
            sensor_config.gps_position_sigma.setZero(); sensor_config.gps_velocity_sigma.setZero();
            sensor_config.baro_pressure_sigma_pa=0.0;
        }
        sensor_config.accel_bias_sigma*=std::max(0.0,bias_scale);
        sensor_config.gyro_bias_sigma*=std::max(0.0,bias_scale);
        dynamics::State x; x.position_i=Eigen::Vector3d(p.planet_radius_m+0.01,0,0); x.velocity_i.setZero(); x.propellant_mass=p.initial_propellant_mass_kg;
        Eigen::Matrix<double,estimator::kErrorN,estimator::kErrorN> P=Eigen::Matrix<double,estimator::kErrorN,estimator::kErrorN>::Zero();
        P.block<3,3>(0,0).diagonal().setConstant(3.0*3.0);
        P.block<3,3>(3,3).diagonal().setConstant(0.5*0.5);
        P.block<3,3>(6,6).diagonal().setConstant(std::pow(0.25*0.5*deg2rad,2));
        P.block<3,3>(9,9).diagonal()=sensor_config.accel_bias_sigma.cwiseAbs2();
        P.block<3,3>(12,12).diagonal()=sensor_config.gyro_bias_sigma.cwiseAbs2();
        dynamics::State initial_estimate=x;
        for(int axis=0;axis<3;++axis) {
            initial_estimate.position_i(axis)+=3.0*w(rng);
            initial_estimate.velocity_i(axis)+=0.5*w(rng);
        }
        if(c.navigation_only) initial_estimate=x;
        estimator::UKF filter(initial_estimate,P); sensors::SensorSuite suite(sensor_config,0xabc000ULL+static_cast<unsigned>(i)); sensors::MeasurementQueue queue; control::AttitudeController controller;
        double apogee=x.position_i.norm()-p.planet_radius_m; dynamics::State impact=x; bool landed=false;
        double nees_sum=0.0, nis_sum=0.0, gps_nis_sum=0.0, baro_nis_sum=0.0;
        int nees_count=0, nis_count=0, gps_nis_count=0, baro_nis_count=0;
        double next_consistency_time=0.0;
        for(int k=0;k<static_cast<int>(c.duration_s/c.dt_s);++k) {
            const double t=k*c.dt_s; const auto previous=x; dynamics::Inputs input{t,0,0};
            if(c.closed_loop) { const double h=x.position_i.norm()-p.planet_radius_m; const auto env=dynamics::environment_at(p,h); const auto phase=control::phase_from(t,h,x.velocity_i.dot(x.position_i.normalized())); auto guidance=filter.estimate(); if(c.truth_guidance) guidance.nominal=x; const auto tvc_cmd=controller.command(guidance,phase,0.5*env.density_kgpm3*x.velocity_i.squaredNorm()); const auto tvc=controller.actuator().update(tvc_cmd,c.dt_s); input.tvc_pitch_rad=tvc(0); input.tvc_yaw_rad=tvc(1); }
            const auto dx=dynamics::derivative(x,p,input);
            if(c.closed_loop || c.navigation_only) {
                for(const auto& m:suite.sample(x,dx,t,c.dt_s,p)) queue.push(m);
                for(const auto& m:queue.release(t)) {
                    filter.update(m,p);
                    if(!std::holds_alternative<sensors::ImuMeasurement>(m) &&
                       filter.last_nis_dof()>0) {
                        const double normalized_nis=filter.last_nis()/std::max(1,filter.last_nis_dof());
                        nis_sum += normalized_nis;
                        ++nis_count;
                        if(std::holds_alternative<sensors::GpsMeasurement>(m)) {
                            gps_nis_sum+=normalized_nis; ++gps_nis_count;
                            if(consistency_samples) traces[i].push_back(
                                {i,std::get<sensors::GpsMeasurement>(m).measurement_time_s,
                                 "gps_nis",filter.last_nis(),filter.last_nis_dof()});
                        } else {
                            baro_nis_sum+=normalized_nis; ++baro_nis_count;
                            if(consistency_samples) traces[i].push_back(
                                {i,std::get<sensors::BaroMeasurement>(m).time_s,
                                 "baro_nis",filter.last_nis(),filter.last_nis_dof()});
                        }
                    }
                }
            }
            x=dynamics::rk4_step(x,p,input,c.dt_s); if(c.closed_loop || c.navigation_only) filter.propagate(c.dt_s,p,input);
            if(c.closed_loop || c.navigation_only) {
                Eigen::Matrix<double,6,1> err;
                err<<x.position_i-filter.estimate().nominal.position_i,
                     x.velocity_i-filter.estimate().nominal.velocity_i;
                const auto cov=filter.estimate().covariance.block<6,6>(0,0);
                const double nees=(err.transpose()*cov.ldlt().solve(err))(0);
                if(std::isfinite(nees)) {
                    nees_sum+=nees; ++nees_count;
                    const double state_time=t+c.dt_s;
                    if(consistency_samples && state_time+0.5*c.dt_s>=next_consistency_time) {
                        traces[i].push_back({i,state_time,"nees",nees,6});
                        next_consistency_time += std::max(c.consistency_interval_s,c.dt_s);
                    }
                }
            }
            const double h0=previous.position_i.norm()-p.planet_radius_m, h1=x.position_i.norm()-p.planet_radius_m; apogee=std::max(apogee,h1);
            if(h1<=0.0) { const double f=h0/(h0-h1); impact=previous; impact.position_i=previous.position_i+f*(x.position_i-previous.position_i); impact.velocity_i=previous.velocity_i+f*(x.velocity_i-previous.velocity_i); landed=true; break; }
        }
        const double downrange=landed?impact.position_i.y():std::numeric_limits<double>::quiet_NaN(); const double crossrange=landed?impact.position_i.z():std::numeric_limits<double>::quiet_NaN();
        out[i]={i,apogee,downrange,crossrange,wind_scale,thrust_scale,dry_mass_scale,aero_scale,bias_scale,misalignment_deg,cg_offset_m,nees_count?nees_sum/nees_count:0.0,nis_count?nis_sum/nis_count:0.0,gps_nis_count?gps_nis_sum/gps_nis_count:0.0,baro_nis_count?baro_nis_sum/baro_nis_count:0.0,landed?1:0};
    }
    if(consistency_samples) {
        consistency_samples->clear();
        std::size_t total=0;
        for(const auto& trace:traces) total+=trace.size();
        consistency_samples->reserve(total);
        for(auto& trace:traces)
            consistency_samples->insert(consistency_samples->end(),trace.begin(),trace.end());
    }
    return out;
}
void write_csv(const std::vector<RunResult>& r,const std::string& path) { std::ofstream f(path); f<<"run,apogee_m,impact_x_m,impact_y_m,wind_scale,thrust_scale,dry_mass_scale,aero_scale,sensor_bias_scale,thrust_misalignment_deg,cg_offset_m,mean_nees,mean_nis,mean_gps_nis,mean_baro_nis,impact_valid\n"; for(auto& x:r)f<<x.run<<","<<x.apogee_m<<","<<x.impact_x_m<<","<<x.impact_y_m<<","<<x.wind_scale<<","<<x.thrust_scale<<","<<x.dry_mass_scale<<","<<x.aero_scale<<","<<x.sensor_bias_scale<<","<<x.thrust_misalignment_deg<<","<<x.cg_offset_m<<","<<x.mean_nees<<","<<x.mean_nis<<","<<x.mean_gps_nis<<","<<x.mean_baro_nis<<","<<x.impact_valid<<"\n"; }
void write_consistency_csv(const std::vector<ConsistencySample>& samples,const std::string& path) {
    std::ofstream f(path);
    f<<"run,time_s,metric,value,dof\n";
    for(const auto& x:samples)
        f<<x.run<<","<<x.time_s<<","<<x.metric<<","<<x.value<<","<<x.dof<<"\n";
}
}
