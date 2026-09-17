#pragma once
#include "estimator/ukf.hpp"

namespace lv::control {
enum class Phase { Launch, GravityTurn, Coast, Descent, Landed };
struct ActuatorConfig { double lag_s{0.08}; double slew_radps{0.25}; double limit_rad{0.08}; };
class TvcrActuator {
public:
    explicit TvcrActuator(ActuatorConfig c={}) : c_(c) {}
    Eigen::Vector2d update(const Eigen::Vector2d& command, double dt);
    const Eigen::Vector2d& state() const { return state_; }
private: ActuatorConfig c_; Eigen::Vector2d state_{Eigen::Vector2d::Zero()};
};
class AttitudeController {
public:
    AttitudeController();
    Eigen::Vector2d command(const estimator::Estimate& estimate, Phase phase, double dynamic_pressure_pa);
    TvcrActuator& actuator() { return actuator_; }
private: Eigen::Matrix<double,2,6> k_launch_, k_turn_, k_coast_; TvcrActuator actuator_;
};
Phase phase_from(double t, double altitude_m, double vertical_speed_mps);
}
