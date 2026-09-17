#include "control/control.hpp"
#include <algorithm>
#include <cmath>
namespace lv::control {
Eigen::Vector2d TvcrActuator::update(const Eigen::Vector2d& cmd, double dt) {
    Eigen::Vector2d target=cmd.cwiseMax(-c_.limit_rad).cwiseMin(c_.limit_rad);
    const double a=dt/(c_.lag_s+dt); Eigen::Vector2d desired=state_+a*(target-state_);
    Eigen::Vector2d d=desired-state_; const double lim=c_.slew_radps*dt;
    for(int i=0;i<2;++i) d(i)=std::clamp(d(i),-lim,lim); state_+=d; return state_;
}
AttitudeController::AttitudeController() : actuator_({0.03,5.0,0.08}) {
    // The vehicle thrust axis is body +X.  Its pointing error therefore lives
    // in the inertial Y/Z components of x_body.cross(x_desired), not X.
    // Rows are TVC pitch/yaw; columns are [attitude error(3), rate(3)].
    k_launch_<<0.0,0.7,0.0, 0.0,0.15,0.0,
               0.0,0.0,0.7, 0.0,0.0,0.15;
    k_turn_=k_launch_; k_coast_=0.5*k_launch_;
}
Eigen::Vector2d AttitudeController::command(const estimator::Estimate& e, Phase phase, double q) {
    const auto& K=(q<100.0)?k_launch_:((phase==Phase::Coast)?k_coast_:k_turn_);
    const Eigen::Vector3d thrust_axis_i=e.nominal.attitude_bi*Eigen::Vector3d::UnitX();
    // With the nozzle arm aft of the CG, this cross-product orientation makes
    // -K*error generate a restoring TVC moment.
    const Eigen::Vector3d pointing_error=thrust_axis_i.cross(Eigen::Vector3d::UnitX());
    Eigen::Matrix<double,6,1> z; z<<pointing_error,e.nominal.angular_rate_b;
    // The dynamics convention uses q_BI and an aft nozzle arm; under that
    // convention the restoring command has the positive feedback sign here.
    return K*z;
}
Phase phase_from(double t,double h,double vz) { if(h<10 && t<5)return Phase::Launch; if(vz>0)return Phase::GravityTurn; if(vz<0)return Phase::Descent; return Phase::Coast; }
}
