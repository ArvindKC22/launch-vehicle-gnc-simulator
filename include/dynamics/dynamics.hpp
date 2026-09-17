#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace lv::dynamics {

using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;

// Quaternion convention: q_BI rotates body-frame vectors into the inertial frame.
// State ordering is [r_I(3), v_I(3), q_BI(w,x,y,z), omega_B(3), propellant_mass].
struct State {
    Vec3 position_i{Vec3::Zero()};
    Vec3 velocity_i{Vec3::Zero()};
    Eigen::Quaterniond attitude_bi{Eigen::Quaterniond::Identity()};
    Vec3 angular_rate_b{Vec3::Zero()};
    double propellant_mass{0.0};
};

struct StateDerivative {
    Vec3 position_i_dot{Vec3::Zero()};
    Vec3 velocity_i_dot{Vec3::Zero()};
    Eigen::Quaterniond attitude_bi_dot{0.0, 0.0, 0.0, 0.0};
    Vec3 angular_rate_b_dot{Vec3::Zero()};
    double propellant_mass_dot{0.0};
};

struct TablePoint { double x; double value; };

struct WindPoint { double altitude_m; Vec3 velocity_i_mps; };

struct AeroTable {
    std::vector<TablePoint> lift_vs_alpha;
    std::vector<TablePoint> drag_vs_mach;
};

struct MotorPoint {
    double time_s;
    double vacuum_thrust_n;
    double mass_flow_kgps;
};

struct VehicleParameters {
    double dry_mass_kg{0.0};
    double initial_propellant_mass_kg{0.0};
    double reference_area_m2{0.0};
    // Vector from CG to the thrust application point, expressed in body axes.
    Vec3 thrust_application_point_b_m{Vec3::Zero()};
    double nozzle_exit_area_m2{0.0};
    double thrust_misalignment_pitch_rad{0.0};
    double thrust_misalignment_yaw_rad{0.0};
    Mat3 inertia_dry_b{Mat3::Identity()};
    Mat3 inertia_full_b{Mat3::Identity()};
    double gravitational_parameter_m3ps2{3.986004418e14};
    double planet_radius_m{6378137.0};
    double atmosphere_max_altitude_m{200000.0};
    AeroTable aero;
    std::vector<MotorPoint> motor;
    std::vector<WindPoint> wind;
};

struct Environment {
    double pressure_pa{0.0};
    double density_kgpm3{0.0};
    double speed_of_sound_mps{340.0};
    Vec3 wind_i_mps{Vec3::Zero()};
};

struct Inputs {
    double time_s{0.0};
    double tvc_pitch_rad{0.0};
    double tvc_yaw_rad{0.0};
};

Environment us_standard_atmosphere(double altitude_m);
Environment environment_at(const VehicleParameters& p, double altitude_m);
double interpolate(const std::vector<TablePoint>& table, double x, double fallback = 0.0);

double mass(const State& x, const VehicleParameters& p);
Mat3 inertia_body(const State& x, const VehicleParameters& p);
StateDerivative derivative(const State& x, const VehicleParameters& p, const Inputs& u);
State rk4_step(const State& x, const VehicleParameters& p, const Inputs& u, double dt_s);

// Specific mechanical energy for the central-gravity ballistic validation.
double specific_orbital_energy(const State& x, const VehicleParameters& p);

}  // namespace lv::dynamics
