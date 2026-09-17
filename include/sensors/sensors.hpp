#pragma once

#include "dynamics/dynamics.hpp"
#include <random>
#include <cstdint>
#include <variant>
#include <queue>

namespace lv::sensors {

struct SensorConfig {
    double imu_rate_hz{200.0};
    double gps_rate_hz{10.0};
    double gps_latency_s{0.20};
    Eigen::Vector3d accel_bias_sigma{Eigen::Vector3d::Constant(0.03)};
    Eigen::Vector3d gyro_bias_sigma{Eigen::Vector3d::Constant(0.002)};
    Eigen::Vector3d accel_noise_sigma{Eigen::Vector3d::Constant(0.02)};
    Eigen::Vector3d gyro_noise_sigma{Eigen::Vector3d::Constant(0.0005)};
    Eigen::Vector3d gps_position_sigma{Eigen::Vector3d::Constant(3.0)};
    Eigen::Vector3d gps_velocity_sigma{Eigen::Vector3d::Constant(0.15)};
    double baro_rate_hz{50.0};
    double baro_lag_s{0.15};
    double baro_pressure_sigma_pa{25.0};
    Eigen::Vector3d accel_scale_error{Eigen::Vector3d::Zero()};
    Eigen::Vector3d gyro_scale_error{Eigen::Vector3d::Zero()};
    Eigen::Vector3d accel_bias_rw_sigma{Eigen::Vector3d::Constant(0.0002)};
    Eigen::Vector3d gyro_bias_rw_sigma{Eigen::Vector3d::Constant(0.00001)};
};

struct ImuMeasurement { double time_s; Eigen::Vector3d specific_force_b; Eigen::Vector3d angular_rate_b; };
struct BaroMeasurement { double time_s; double altitude_m; };
struct GpsMeasurement { double measurement_time_s; double arrival_time_s; Eigen::Vector3d position_i; Eigen::Vector3d velocity_i; };
using Measurement = std::variant<ImuMeasurement, BaroMeasurement, GpsMeasurement>;

struct QueuedMeasurement {
    double arrival_time_s;
    Measurement measurement;
    bool operator>(const QueuedMeasurement& other) const { return arrival_time_s > other.arrival_time_s; }
};

class MeasurementQueue {
public:
    void push(const Measurement& measurement);
    std::vector<Measurement> release(double current_time_s);
    std::size_t size() const { return queue_.size(); }
private:
    std::priority_queue<QueuedMeasurement, std::vector<QueuedMeasurement>, std::greater<QueuedMeasurement>> queue_;
};

class SensorSuite {
public:
    SensorSuite(SensorConfig config, std::uint64_t seed = 1);
    std::vector<Measurement> sample(const dynamics::State& truth, const dynamics::StateDerivative& truth_dot,
                                    double time_s, double dt_s, const dynamics::VehicleParameters& vehicle);
    const Eigen::Vector3d& accel_bias() const { return accel_bias_; }
    const Eigen::Vector3d& gyro_bias() const { return gyro_bias_; }
private:
    double next_imu_{0.0}, next_gps_{0.0}, next_baro_{0.0};
    SensorConfig c_;
    std::mt19937_64 rng_;
    Eigen::Vector3d accel_bias_{Eigen::Vector3d::Zero()}, gyro_bias_{Eigen::Vector3d::Zero()};
    double baro_state_{0.0};
    Eigen::Vector3d normal3(const Eigen::Vector3d& sigma);
};

} // namespace lv::sensors
