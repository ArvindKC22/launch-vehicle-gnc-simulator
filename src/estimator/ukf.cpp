#include "estimator/ukf.hpp"
#include <Eigen/Dense>

namespace lv::estimator {

Eigen::Vector3d mrp_shadow_set(const Eigen::Vector3d& sigma) {
    const double s2=sigma.squaredNorm();
    return (s2>1.0) ? -sigma/s2 : sigma;
}

namespace {
Eigen::Vector3d mrp_from_quat(const Eigen::Quaterniond& q) {
    Eigen::Quaterniond n=q; n.normalize();
    const double den=1.0+n.w();
    if (std::abs(den)<1e-10) return Eigen::Vector3d::Zero();
    Eigen::Vector3d sigma=n.vec()/den; return mrp_shadow_set(sigma);
}
Eigen::Quaterniond quat_from_mrp(const Eigen::Vector3d& s) {
    const double s2=s.squaredNorm();
    return Eigen::Quaterniond((1-s2)/(1+s2), 2*s.x()/(1+s2), 2*s.y()/(1+s2), 2*s.z()/(1+s2));
}
void inject_estimate(Estimate& e, const Eigen::VectorXd& d) {
    e.nominal.position_i += d.segment<3>(0);
    e.nominal.velocity_i += d.segment<3>(3);
    const Eigen::Quaterniond dq=quat_from_mrp(mrp_shadow_set(d.segment<3>(6)));
    e.nominal.attitude_bi=(e.nominal.attitude_bi*dq).normalized();
    e.accel_bias += d.segment<3>(9);
    e.gyro_bias += d.segment<3>(12);
}

Eigen::Matrix<double,kErrorN,kErrorN> process_noise(double dt, bool powered) {
    Eigen::Matrix<double,kErrorN,kErrorN> q=
        Eigen::Matrix<double,kErrorN,kErrorN>::Zero();
    // White acceleration model.  Position/velocity terms are the exact
    // discrete covariance for acceleration held over one sample.
    // IMU propagation intentionally does not assume the exact thrust,
    // aerodynamic, or mass model.  The powered value absorbs launch-vehicle
    // model error; the coast value covers residual accelerometer/attitude
    // uncertainty without making the filter needlessly diffuse.
    const double accel_model_sigma=powered ? 5.0 : 160.0;
    const double a2=accel_model_sigma*accel_model_sigma;
    for(int axis=0;axis<3;++axis) {
        q(axis,axis)=0.25*dt*dt*dt*dt*a2;
        q(axis,axis+3)=q(axis+3,axis)=0.5*dt*dt*dt*a2;
        q(axis+3,axis+3)=dt*dt*a2;
    }
    const double gyro_model_sigma=powered ? 0.50 : 0.01;
    q.block<3,3>(6,6).diagonal().setConstant(
        std::pow(0.25*gyro_model_sigma*dt,2));
    q.block<3,3>(9,9).diagonal().setConstant(0.0002*0.0002*dt);
    q.block<3,3>(12,12).diagonal().setConstant(0.00001*0.00001*dt);
    return q;
}

Eigen::Vector3d gravity_at(const dynamics::State& x,const dynamics::VehicleParameters& v) {
    const double r=x.position_i.norm();
    return -v.gravitational_parameter_m3ps2*x.position_i/std::pow(r,3);
}
dynamics::StateDerivative imu_derivative(const dynamics::State& x,const dynamics::VehicleParameters& v,
                                         const Eigen::Vector3d& f_b,const Eigen::Vector3d& w_b,
                                         const Eigen::Vector3d& accel_bias,const Eigen::Vector3d& gyro_bias) {
    dynamics::StateDerivative d;
    d.position_i_dot=x.velocity_i;
    d.velocity_i_dot=gravity_at(x,v)+x.attitude_bi.toRotationMatrix()*(f_b-accel_bias);
    const Eigen::Quaterniond omega_q(0.0,(w_b-gyro_bias).x(),(w_b-gyro_bias).y(),(w_b-gyro_bias).z());
    d.attitude_bi_dot.coeffs()=(x.attitude_bi*omega_q).coeffs()*0.5;
    d.angular_rate_b_dot.setZero();
    d.propellant_mass_dot=0.0;
    return d;
}
dynamics::State imu_rk4_step(const dynamics::State& x,const dynamics::VehicleParameters& v,
                              const Eigen::Vector3d& f_b,const Eigen::Vector3d& w_b,
                              const Eigen::Vector3d& accel_bias,const Eigen::Vector3d& gyro_bias,double dt) {
    const auto k1=imu_derivative(x,v,f_b,w_b,accel_bias,gyro_bias);
    auto add=[&](const dynamics::State& s,const dynamics::StateDerivative& d,double h) {
        auto y=s; y.position_i+=h*d.position_i_dot; y.velocity_i+=h*d.velocity_i_dot;
        y.attitude_bi.coeffs()+=h*d.attitude_bi_dot.coeffs(); y.angular_rate_b+=h*d.angular_rate_b_dot;
        y.propellant_mass+=h*d.propellant_mass_dot; y.attitude_bi.normalize(); return y;
    };
    const auto k2=imu_derivative(add(x,k1,0.5*dt),v,f_b,w_b,accel_bias,gyro_bias);
    const auto k3=imu_derivative(add(x,k2,0.5*dt),v,f_b,w_b,accel_bias,gyro_bias);
    const auto k4=imu_derivative(add(x,k3,dt),v,f_b,w_b,accel_bias,gyro_bias);
    auto y=x; y.position_i+=(dt/6.0)*(k1.position_i_dot+2*k2.position_i_dot+2*k3.position_i_dot+k4.position_i_dot);
    y.velocity_i+=(dt/6.0)*(k1.velocity_i_dot+2*k2.velocity_i_dot+2*k3.velocity_i_dot+k4.velocity_i_dot);
    y.attitude_bi.coeffs()+=(dt/6.0)*(k1.attitude_bi_dot.coeffs()+2*k2.attitude_bi_dot.coeffs()+2*k3.attitude_bi_dot.coeffs()+k4.attitude_bi_dot.coeffs());
    y.attitude_bi.normalize();
    y.angular_rate_b=w_b-gyro_bias;
    return y;
}

void enforce_covariance(Eigen::Matrix<double,kErrorN,kErrorN>& p) {
    p=0.5*(p+p.transpose());
    p.diagonal().array() += 1e-12;
}
}
LinearFilterResult linear_ukf_step(const Eigen::VectorXd& x, const Eigen::MatrixXd& P,
                                   const Eigen::MatrixXd& F, const Eigen::MatrixXd& Q,
                                   const Eigen::MatrixXd& H, const Eigen::MatrixXd& R,
                                   const Eigen::VectorXd& z) {
    // alpha=1 avoids unnecessary cancellation in this linear reference case;
    // the nonlinear flight UKF uses its separately tuned alpha value.
    const int n=x.size(); const double alpha=1.0, beta=2.0, kappa=0.0;
    const double lambda=alpha*alpha*(n+kappa)-n; const int ns=2*n+1;
    Eigen::VectorXd wm=Eigen::VectorXd::Zero(ns),wc=Eigen::VectorXd::Zero(ns);
    wm(0)=lambda/(n+lambda); wc(0)=wm(0)+1-alpha*alpha+beta;
    for(int i=1;i<ns;++i) { wm(i)=1.0/(2*(n+lambda)); wc(i)=wm(i); }
    // Factor P first, then apply the sigma-point scale explicitly. This keeps
    // the factorization numerically transparent and avoids scaled-expression
    // lifetime/triangular-view surprises across Eigen versions.
    Eigen::MatrixXd L=P.llt().matrixL().toDenseMatrix();
    L.triangularView<Eigen::StrictlyUpper>().setZero();
    L*=std::sqrt(n+lambda); Eigen::MatrixXd X(n,ns); X.col(0)=x;
    for(int i=0;i<n;++i){X.col(i+1)=x+L.col(i); X.col(i+1+n)=x-L.col(i);}
    Eigen::MatrixXd Y=F*X; Eigen::VectorXd xp=Eigen::VectorXd::Zero(n); for(int i=0;i<ns;++i)xp+=wm(i)*Y.col(i);
    Eigen::MatrixXd Pp=Q; for(int i=0;i<ns;++i){auto d=Y.col(i)-xp; Pp+=wc(i)*(d*d.transpose());}
    Eigen::MatrixXd Z=H*Y; const int m=z.size(); Eigen::VectorXd zp=Eigen::VectorXd::Zero(m); for(int i=0;i<ns;++i)zp+=wm(i)*Z.col(i);
    // Y contains the deterministic F*x sigma points. Process noise was added
    // to Pp, so its measurement-space contribution must also enter S.
    Eigen::MatrixXd S=H*Q*H.transpose()+R, Pxz=Q*H.transpose(); for(int i=0;i<ns;++i){auto dz=Z.col(i)-zp; auto dx=Y.col(i)-xp; S+=wc(i)*(dz*dz.transpose()); Pxz+=wc(i)*(dx*dz.transpose());}
    Eigen::VectorXd innov=z-zp; Eigen::MatrixXd K=Pxz*S.ldlt().solve(Eigen::MatrixXd::Identity(m,m)); return {xp+K*innov,Pp-K*S*K.transpose()};
}
UKF::UKF(const dynamics::State& initial, const Eigen::Matrix<double,kErrorN,kErrorN>& covariance) {
    e_.nominal=initial; e_.nominal.attitude_bi.normalize(); e_.covariance=covariance;
    state_history_.push_back({0.0,e_.nominal.position_i,e_.nominal.velocity_i});
}
void UKF::inject(const Eigen::VectorXd& d) {
    inject_estimate(e_, d);
}
void UKF::propagate(double dt, const dynamics::VehicleParameters& v, const dynamics::Inputs& u) {
    last_input_=u;
    const int n=kErrorN; const double lambda=alpha_*alpha_*(n+kappa_)-n; const int m=2*n+1;
    enforce_covariance(e_.covariance);
    Eigen::LLT<Eigen::Matrix<double,kErrorN,kErrorN>> llt(e_.covariance);
    if(llt.info()!=Eigen::Success) {
        e_.covariance.diagonal().array() += 1e-8;
        llt.compute(e_.covariance);
    }
    if(llt.info()!=Eigen::Success) throw std::runtime_error("UKF covariance is not positive definite");
    Eigen::MatrixXd L=llt.matrixL().toDenseMatrix(); L.triangularView<Eigen::StrictlyUpper>().setZero(); L*=std::sqrt(n+lambda);
    std::vector<Estimate> sig(m, e_); std::vector<double> wm(m), wc(m); wm[0]=lambda/(n+lambda); wc[0]=wm[0]+1-alpha_*alpha_+beta_;
    for(int i=1;i<m;++i) { wm[i]=1.0/(2*(n+lambda)); wc[i]=wm[i]; }
    for(int i=0;i<n;++i) { sig[1+i]=e_; inject_estimate(sig[1+i], L.col(i)); }
    for(int i=0;i<n;++i) { sig[1+n+i]=e_; inject_estimate(sig[1+n+i], -L.col(i)); }
    for(int i=0;i<m;++i) {
        if(have_imu_measurement_)
            sig[i].nominal=imu_rk4_step(sig[i].nominal,v,last_specific_force_b_,
                                        last_angular_rate_b_,sig[i].accel_bias,
                                        sig[i].gyro_bias,dt);
        else
            sig[i].nominal=dynamics::rk4_step(sig[i].nominal,v,u,dt);
        // Propellant depletion is not observable from the IMU kinematics, but
        // it is a known commanded-time state in this vehicle model.  Preserve
        // that state in the IMU-driven branch so powered/coast noise tuning
        // and inertia scheduling switch at the actual burnout epoch.
        if(have_imu_measurement_) {
            const auto model_rate=dynamics::derivative(sig[i].nominal,v,u).propellant_mass_dot;
            sig[i].nominal.propellant_mass=std::clamp(
                sig[i].nominal.propellant_mass+dt*model_rate,0.0,v.initial_propellant_mass_kg);
        }
    }
    Estimate mean=sig[0];
    mean.nominal.position_i.setZero(); mean.nominal.velocity_i.setZero();
    mean.nominal.angular_rate_b.setZero();
    mean.accel_bias.setZero(); mean.gyro_bias.setZero();
    for(int i=0;i<m;++i) {
        mean.nominal.position_i += wm[i]*sig[i].nominal.position_i;
        mean.nominal.velocity_i += wm[i]*sig[i].nominal.velocity_i;
        mean.nominal.angular_rate_b += wm[i]*sig[i].nominal.angular_rate_b;
        mean.accel_bias += wm[i]*sig[i].accel_bias;
        mean.gyro_bias += wm[i]*sig[i].gyro_bias;
    }
    // Compute the quaternion mean on the local MRP tangent space.  This
    // avoids arithmetic quaternion averaging and remains in the shadow set.
    mean.nominal.attitude_bi=sig[0].nominal.attitude_bi.normalized();
    for(int iteration=0;iteration<8;++iteration) {
        Eigen::Vector3d correction=Eigen::Vector3d::Zero();
        for(int i=0;i<m;++i)
            correction += wm[i]*mrp_from_quat(
                mean.nominal.attitude_bi.conjugate()*sig[i].nominal.attitude_bi);
        if(correction.norm()<1e-13) break;
        mean.nominal.attitude_bi=(mean.nominal.attitude_bi*quat_from_mrp(correction)).normalized();
    }
    mean.covariance.setZero();
    for(int i=0;i<m;++i) {
        Eigen::Matrix<double,kErrorN,1> d;
        d<<sig[i].nominal.position_i-mean.nominal.position_i,
           sig[i].nominal.velocity_i-mean.nominal.velocity_i,
           mrp_from_quat(mean.nominal.attitude_bi.conjugate()*sig[i].nominal.attitude_bi),
           sig[i].accel_bias-mean.accel_bias,
           sig[i].gyro_bias-mean.gyro_bias;
        mean.covariance += wc[i]*(d*d.transpose());
    }
    const bool powered=!v.motor.empty() && e_.nominal.propellant_mass>0.0 &&
        u.time_s<v.motor.back().time_s;
    mean.covariance += process_noise(dt,powered);
    enforce_covariance(mean.covariance);
    e_=mean;
    current_time_s_ += dt;
    state_history_.push_back(
        {current_time_s_,e_.nominal.position_i,e_.nominal.velocity_i});
    while(state_history_.size()>2 && state_history_[1].time_s<current_time_s_-2.0)
        state_history_.pop_front();
}
std::pair<Eigen::Vector3d,Eigen::Vector3d> UKF::kinematics_at(double time_s) const {
    if(state_history_.empty() || time_s>=state_history_.back().time_s)
        return {e_.nominal.position_i,e_.nominal.velocity_i};
    if(time_s<=state_history_.front().time_s)
        return {state_history_.front().position_i,state_history_.front().velocity_i};
    const auto upper=std::lower_bound(
        state_history_.begin(),state_history_.end(),time_s,
        [](const TimedKinematics& state,double time){ return state.time_s<time; });
    const auto lower=std::prev(upper);
    const double fraction=(time_s-lower->time_s)/(upper->time_s-lower->time_s);
    return {lower->position_i+fraction*(upper->position_i-lower->position_i),
            lower->velocity_i+fraction*(upper->velocity_i-lower->velocity_i)};
}
void UKF::update(const sensors::Measurement& z, const dynamics::VehicleParameters& v) {
    Eigen::VectorXd innovation; Eigen::MatrixXd H; Eigen::MatrixXd R;
    const auto* gps_measurement=std::get_if<sensors::GpsMeasurement>(&z);
    if(const auto* imu=std::get_if<sensors::ImuMeasurement>(&z)) {
        innovation=imu->angular_rate_b-(e_.nominal.angular_rate_b+e_.gyro_bias);
        H=Eigen::MatrixXd::Zero(3,kErrorN);
        H.block<3,3>(0,12).setIdentity();
        R=Eigen::Matrix3d::Identity()*0.0005*0.0005;
        last_imu_time_s_=imu->time_s;
        last_specific_force_b_=imu->specific_force_b;
        last_angular_rate_b_=imu->angular_rate_b;
        have_imu_measurement_=true;
    } else if(const auto* g=std::get_if<sensors::GpsMeasurement>(&z)) {
        // Evaluate delayed GPS against the predicted acquisition-epoch state.
        // History remains valid through thrust ramps and high-jerk intervals.
        const auto delayed=kinematics_at(g->measurement_time_s);
        innovation.resize(6);
        innovation<<g->position_i-delayed.first,g->velocity_i-delayed.second;
        H=Eigen::MatrixXd::Zero(6,kErrorN);
        H.block<3,3>(0,0).setIdentity();
        H.block<3,3>(3,3).setIdentity();
        R=Eigen::MatrixXd::Zero(6,6);
        R.block<3,3>(0,0).diagonal().setConstant(9.0);
        R.block<3,3>(3,3).diagonal().setConstant(0.0225);
    } else {
        const auto* b=std::get_if<sensors::BaroMeasurement>(&z);
        const double baro_dt=last_baro_time_s_<0.0 ? 1.0/50.0 :
            std::max(1e-6,b->time_s-last_baro_time_s_);
        constexpr double baro_lag_s=0.15;
        const double alpha=baro_dt/(baro_lag_s+baro_dt);
        const Eigen::Vector3d radial=e_.nominal.position_i.normalized();
        const double geometric_altitude=e_.nominal.position_i.norm()-v.planet_radius_m;
        baro_predicted_altitude_m_ += alpha*(geometric_altitude-baro_predicted_altitude_m_);
        last_baro_time_s_=b->time_s;
        const auto environment=dynamics::environment_at(v,geometric_altitude);
        // Once pressure noise is a substantial fraction of ambient pressure,
        // converting it to altitude is strongly non-Gaussian.  Keep advancing
        // the lag state, but reject that measurement and rely on GPS/IMU.
        if(environment.pressure_pa<500.0) {
            last_nis_=0.0; last_nis_dof_=0; return;
        }
        innovation.resize(1);
        innovation(0)=b->altitude_m-baro_predicted_altitude_m_;
        H=Eigen::MatrixXd::Zero(1,kErrorN);
        H.block<1,3>(0,0)=alpha*radial.transpose();
        // Hydrostatic conversion: sigma_h ~= (R*T/g)*sigma_p/p.  A 6.3 km
        // scale height is representative of the modeled lower atmosphere.
        const double pressure_altitude_sigma=
            6300.0*25.0/environment.pressure_pa;
        // The sensor's first-order lag also low-pass filters its white noise.
        // For y[k]=(1-alpha)y[k-1]+alpha*n[k], the steady-state variance is
        // alpha/(2-alpha) times the input variance.
        const double filtered_noise_variance=
            alpha/(2.0-alpha)*pressure_altitude_sigma*pressure_altitude_sigma;
        constexpr double residual_model_variance=0.5*0.5;
        R=Eigen::MatrixXd::Constant(1,1,
            filtered_noise_variance+residual_model_variance);
    }
    const Eigen::MatrixXd S=H*e_.covariance*H.transpose()+R;
    const Eigen::LDLT<Eigen::MatrixXd> solver(S);
    if(solver.info()!=Eigen::Success) throw std::runtime_error("UKF innovation covariance solve failed");
    const Eigen::MatrixXd K=e_.covariance*H.transpose()*solver.solve(
        Eigen::MatrixXd::Identity(S.rows(),S.cols()));
    last_nis_=(innovation.transpose()*solver.solve(innovation))(0);
    last_nis_dof_=static_cast<int>(innovation.size());
    const Eigen::Vector3d position_before=e_.nominal.position_i;
    const Eigen::Vector3d velocity_before=e_.nominal.velocity_i;
    Eigen::VectorXd correction=K*innovation;
    if(gps_measurement) {
        const double age=std::max(0.0,current_time_s_-gps_measurement->measurement_time_s);
        correction.segment<3>(0)+=age*correction.segment<3>(3);
    }
    inject(correction);
    if(gps_measurement) {
        const Eigen::Vector3d delta_position=e_.nominal.position_i-position_before;
        const Eigen::Vector3d delta_velocity=e_.nominal.velocity_i-velocity_before;
        const double measurement_age=current_time_s_-gps_measurement->measurement_time_s;
        for(auto& history:state_history_) {
            if(history.time_s+1e-9<gps_measurement->measurement_time_s) continue;
            const double elapsed=history.time_s-gps_measurement->measurement_time_s;
            history.position_i+=delta_position+(elapsed-measurement_age)*delta_velocity;
            history.velocity_i+=delta_velocity;
        }
    }
    const Eigen::Matrix<double,kErrorN,kErrorN> identity=
        Eigen::Matrix<double,kErrorN,kErrorN>::Identity();
    const auto a=identity-K*H;
    e_.covariance=a*e_.covariance*a.transpose()+K*R*K.transpose();
    enforce_covariance(e_.covariance);
    if(!state_history_.empty() &&
       std::abs(state_history_.back().time_s-current_time_s_)<1e-9) {
        state_history_.back().position_i=e_.nominal.position_i;
        state_history_.back().velocity_i=e_.nominal.velocity_i;
    }
}

EKF::EKF(const dynamics::State& initial, const Eigen::Matrix<double,kErrorN,kErrorN>& covariance) {
    e_.nominal=initial; e_.nominal.attitude_bi.normalize(); e_.covariance=covariance;
}
void EKF::propagate(double dt, const dynamics::VehicleParameters& v, const dynamics::Inputs& u) {
    e_.nominal=dynamics::rk4_step(e_.nominal,v,u,dt);
    e_.covariance += dt*1e-6*Eigen::Matrix<double,kErrorN,kErrorN>::Identity();
}
void EKF::update(const sensors::Measurement& z, const dynamics::VehicleParameters& v) {
    if(std::holds_alternative<sensors::ImuMeasurement>(z)) return;
    Eigen::MatrixXd H; Eigen::VectorXd innovation; Eigen::MatrixXd R;
    if(auto* g=std::get_if<sensors::GpsMeasurement>(&z)) {
        H=Eigen::MatrixXd::Zero(6,kErrorN); H.block(0,0,3,3).setIdentity(); H.block(3,3,3,3).setIdentity(); innovation.resize(6); innovation<<g->position_i-e_.nominal.position_i,g->velocity_i-e_.nominal.velocity_i; R=Eigen::MatrixXd::Identity(6,6); R.diagonal()<<9,9,9,.0225,.0225,.0225;
    } else {
        auto* b=std::get_if<sensors::BaroMeasurement>(&z); H=Eigen::MatrixXd::Zero(1,kErrorN); H.block<1,3>(0,0)=e_.nominal.position_i.normalized().transpose(); innovation.resize(1); innovation(0)=b->altitude_m-(e_.nominal.position_i.norm()-v.planet_radius_m); R=Eigen::MatrixXd::Constant(1,1,625.0);
    }
    const Eigen::MatrixXd S=H*e_.covariance*H.transpose()+R; const Eigen::MatrixXd K=e_.covariance*H.transpose()*S.ldlt().solve(Eigen::MatrixXd::Identity(S.rows(),S.cols())); last_nis_=(innovation.transpose()*S.ldlt().solve(innovation))(0); e_.nominal.position_i+= (K*innovation).segment<3>(0); e_.nominal.velocity_i+=(K*innovation).segment<3>(3); e_.covariance=(Eigen::Matrix<double,kErrorN,kErrorN>::Identity()-K*H)*e_.covariance;
}
} // namespace lv::estimator
