#pragma once

#include "dynamics/dynamics.hpp"
#include "sensors/sensors.hpp"
#include <Eigen/Cholesky>
#include <deque>
#include <functional>
#include <utility>

namespace lv::estimator {

// Error state: [position error, velocity error, attitude MRP error,
// accelerometer bias, gyro bias]. Quaternion is stored separately in nominal_.
constexpr int kErrorN = 15;
struct Estimate {
    dynamics::State nominal;
    Eigen::Vector3d accel_bias{Eigen::Vector3d::Zero()};
    Eigen::Vector3d gyro_bias{Eigen::Vector3d::Zero()};
    Eigen::Matrix<double,kErrorN,kErrorN> covariance;
};

struct LinearFilterResult { Eigen::VectorXd state; Eigen::MatrixXd covariance; };
Eigen::Vector3d mrp_shadow_set(const Eigen::Vector3d& sigma);
// Reference implementation used only for validation of the unscented weights.
// For linear models it must agree with the analytic Kalman filter.
LinearFilterResult linear_ukf_step(const Eigen::VectorXd& state, const Eigen::MatrixXd& covariance,
                                   const Eigen::MatrixXd& F, const Eigen::MatrixXd& Q,
                                   const Eigen::MatrixXd& H, const Eigen::MatrixXd& R,
                                   const Eigen::VectorXd& measurement);

class UKF {
public:
    UKF(const dynamics::State& initial, const Eigen::Matrix<double,kErrorN,kErrorN>& covariance);
    void propagate(double dt, const dynamics::VehicleParameters& vehicle, const dynamics::Inputs& input);
    void update(const sensors::Measurement& measurement, const dynamics::VehicleParameters& vehicle);
    const Estimate& estimate() const { return e_; }
    double last_nis() const { return last_nis_; }
    int last_nis_dof() const { return last_nis_dof_; }
private:
    Estimate e_{}; double last_nis_{0.0}; int last_nis_dof_{0};
    static constexpr double alpha_=0.35, beta_=2.0, kappa_=0.0;
    void inject(const Eigen::VectorXd& delta);
    double last_imu_time_s_{-1.0};
    double current_time_s_{0.0};
    dynamics::Inputs last_input_{};
    bool have_imu_measurement_{false};
    Eigen::Vector3d last_specific_force_b_{Eigen::Vector3d::Zero()};
    Eigen::Vector3d last_angular_rate_b_{Eigen::Vector3d::Zero()};
    struct TimedKinematics {
        double time_s;
        Eigen::Vector3d position_i;
        Eigen::Vector3d velocity_i;
    };
    std::deque<TimedKinematics> state_history_;
    std::pair<Eigen::Vector3d,Eigen::Vector3d> kinematics_at(double time_s) const;
    double baro_predicted_altitude_m_{0.0};
    double last_baro_time_s_{-1.0};
};

// Separate additive-error EKF with the same state and measurement conventions.
class EKF {
public:
    EKF(const dynamics::State& initial, const Eigen::Matrix<double,kErrorN,kErrorN>& covariance);
    void propagate(double dt, const dynamics::VehicleParameters& vehicle, const dynamics::Inputs& input);
    void update(const sensors::Measurement& measurement, const dynamics::VehicleParameters& vehicle);
    const Estimate& estimate() const { return e_; }
    double last_nis() const { return last_nis_; }
private:
    Estimate e_{}; double last_nis_{0.0};
};

} // namespace lv::estimator
