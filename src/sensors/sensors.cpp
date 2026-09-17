#include "sensors/sensors.hpp"
#include <cmath>

namespace lv::sensors {
namespace {
double pressure_to_standard_altitude(double pressure_pa) {
    constexpr double gas_constant=287.05287;
    constexpr double gravity=9.80665;
    constexpr double lapse_rate=0.0065;
    constexpr double stratosphere_temperature=216.65;
    const double pressure=std::max(pressure_pa,1.0);
    if(pressure>=22632.06) {
        const double temperature=288.15*std::pow(
            pressure/101325.0,gas_constant*lapse_rate/gravity);
        return (288.15-temperature)/lapse_rate;
    }
    if(pressure>=5474.88)
        return 11000.0-gas_constant*stratosphere_temperature/gravity*
            std::log(pressure/22632.06);
    return 20000.0-gas_constant*stratosphere_temperature/gravity*
        std::log(pressure/5474.88);
}
}
SensorSuite::SensorSuite(SensorConfig config, std::uint64_t seed) : c_(config), rng_(seed) {
    // Turn-on biases are one draw per simulated instrument and remain part of
    // the in-run random-walk state thereafter.
    accel_bias_=normal3(c_.accel_bias_sigma);
    gyro_bias_=normal3(c_.gyro_bias_sigma);
}
Eigen::Vector3d SensorSuite::normal3(const Eigen::Vector3d& s) {
    std::normal_distribution<double> n(0.0, 1.0);
    return Eigen::Vector3d(s.x()*n(rng_), s.y()*n(rng_), s.z()*n(rng_));
}
std::vector<Measurement> SensorSuite::sample(const dynamics::State& x, const dynamics::StateDerivative& dx,
                                              double t, double dt, const dynamics::VehicleParameters& p) {
    std::vector<Measurement> out;
    const double h = x.position_i.norm() - p.planet_radius_m;
    const auto env = dynamics::environment_at(p, h);
    const Eigen::Matrix3d R = x.attitude_bi.toRotationMatrix();
    accel_bias_ += normal3(c_.accel_bias_rw_sigma) * std::sqrt(dt);
    gyro_bias_ += normal3(c_.gyro_bias_rw_sigma) * std::sqrt(dt);
    while (t + 1e-12 >= next_imu_) {
        // Specific force excludes gravity; recover it from inertial acceleration.
        const Eigen::Vector3d gravity = -p.gravitational_parameter_m3ps2*x.position_i/std::pow(x.position_i.norm(),3);
        Eigen::Vector3d f_b = R.transpose() * (dx.velocity_i_dot - gravity);
        f_b = (Eigen::Vector3d::Ones()+c_.accel_scale_error).cwiseProduct(f_b + accel_bias_) + normal3(c_.accel_noise_sigma);
        Eigen::Vector3d w = (Eigen::Vector3d::Ones()+c_.gyro_scale_error).cwiseProduct(x.angular_rate_b + gyro_bias_) + normal3(c_.gyro_noise_sigma);
        out.emplace_back(ImuMeasurement{next_imu_, f_b, w});
        next_imu_ += 1.0/c_.imu_rate_hz;
    }
    while (t + 1e-12 >= next_baro_) {
        const double pressure = env.pressure_pa + normal3(Eigen::Vector3d::Constant(c_.baro_pressure_sigma_pa)).x();
        const double measured_h=pressure_to_standard_altitude(pressure);
        const double baro_dt=1.0/c_.baro_rate_hz;
        const double a = baro_dt/(c_.baro_lag_s+baro_dt);
        baro_state_ += a*(measured_h-baro_state_);
        out.emplace_back(BaroMeasurement{next_baro_, baro_state_});
        next_baro_ += 1.0/c_.baro_rate_hz;
    }
    while (t + 1e-12 >= next_gps_) {
        out.emplace_back(GpsMeasurement{next_gps_, next_gps_+c_.gps_latency_s,
            x.position_i + normal3(c_.gps_position_sigma), x.velocity_i + normal3(c_.gps_velocity_sigma)});
        next_gps_ += 1.0/c_.gps_rate_hz;
    }
    return out;
}
} // namespace lv::sensors
