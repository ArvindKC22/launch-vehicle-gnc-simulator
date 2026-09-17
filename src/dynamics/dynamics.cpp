#include "dynamics/dynamics.hpp"

#include <cassert>

namespace lv::dynamics {
namespace {

double clamp01(double x) { return std::clamp(x, 0.0, 1.0); }

double interp(const std::vector<TablePoint>& table, double x, double zero) {
    if (table.empty()) return zero;
    if (x <= table.front().x) return table.front().value;
    if (x >= table.back().x) return table.back().value;
    for (std::size_t i = 1; i < table.size(); ++i) {
        if (x <= table[i].x) {
            const double f = (x - table[i - 1].x) / (table[i].x - table[i - 1].x);
            return table[i - 1].value + f * (table[i].value - table[i - 1].value);
        }
    }
    return table.back().value;
}

Vec3 wind_at(const std::vector<WindPoint>& table, double altitude_m) {
    if (table.empty()) return Vec3::Zero();
    if (altitude_m <= table.front().altitude_m) return table.front().velocity_i_mps;
    if (altitude_m >= table.back().altitude_m) return table.back().velocity_i_mps;
    for (std::size_t i = 1; i < table.size(); ++i) {
        if (altitude_m <= table[i].altitude_m) {
            const double f = (altitude_m - table[i - 1].altitude_m) /
                             (table[i].altitude_m - table[i - 1].altitude_m);
            return (1.0 - f) * table[i - 1].velocity_i_mps + f * table[i].velocity_i_mps;
        }
    }
    return table.back().velocity_i_mps;
}

double motor_value(const std::vector<MotorPoint>& table, double time_s, bool thrust) {
    if (table.empty()) return 0.0;
    const auto value = [thrust](const MotorPoint& m) { return thrust ? m.vacuum_thrust_n : m.mass_flow_kgps; };
    if (time_s <= table.front().time_s) return value(table.front());
    if (time_s >= table.back().time_s) return value(table.back());
    for (std::size_t i=1; i<table.size(); ++i) if (time_s <= table[i].time_s) {
        const double f=(time_s-table[i-1].time_s)/(table[i].time_s-table[i-1].time_s);
        return value(table[i-1]) + f*(value(table[i])-value(table[i-1]));
    }
    return value(table.back());
}

void normalize(State& x) {
    const double n = x.attitude_bi.norm();
    if (!(n > std::numeric_limits<double>::epsilon()))
        throw std::domain_error("attitude quaternion has zero norm");
    x.attitude_bi.coeffs() /= n;
}

State add_scaled(const State& x, const StateDerivative& k, double h) {
    State y = x;
    y.position_i += h * k.position_i_dot;
    y.velocity_i += h * k.velocity_i_dot;
    y.attitude_bi.coeffs() += h * k.attitude_bi_dot.coeffs();
    y.angular_rate_b += h * k.angular_rate_b_dot;
    y.propellant_mass += h * k.propellant_mass_dot;
    return y;
}

}  // namespace

double interpolate(const std::vector<TablePoint>& table, double x, double fallback) {
    return interp(table, x, fallback);
}

Environment us_standard_atmosphere(double h) {
    // 1976 standard atmosphere, troposphere + lower stratosphere. Above 20 km,
    // density is held at the 20 km value until the higher-fidelity atmosphere module.
    constexpr double R = 287.05287;
    constexpr double g0 = 9.80665;
    h = std::max(0.0, h);
    double T, p;
    if (h <= 11000.0) {
        T = 288.15 - 0.0065 * h;
        p = 101325.0 * std::pow(T / 288.15, g0 / (R * 0.0065));
    } else if (h <= 20000.0) {
        T = 216.65;
        p = 22632.06 * std::exp(-g0 * (h - 11000.0) / (R * T));
    } else {
        T = 216.65;
        p = 5474.88 * std::exp(-g0 * (h - 20000.0) / (R * T));
    }
    Environment e;
    e.pressure_pa = p;
    e.density_kgpm3 = p / (R * T);
    e.speed_of_sound_mps = std::sqrt(1.4 * R * T);
    return e;
}

Environment environment_at(const VehicleParameters& p, double altitude_m) {
    Environment e = us_standard_atmosphere(altitude_m);
    e.wind_i_mps = wind_at(p.wind, altitude_m);
    return e;
}

double mass(const State& x, const VehicleParameters& p) {
    return p.dry_mass_kg + std::max(0.0, x.propellant_mass);
}

Mat3 inertia_body(const State& x, const VehicleParameters& p) {
    const double denom = std::max(p.initial_propellant_mass_kg, 1e-12);
    const double f = clamp01(x.propellant_mass / denom);
    return p.inertia_dry_b + f * (p.inertia_full_b - p.inertia_dry_b);
}

StateDerivative derivative(const State& x, const VehicleParameters& p, const Inputs& u) {
    State xn = x;
    normalize(xn);
    const Mat3 r_bi = xn.attitude_bi.toRotationMatrix();
    const double r_norm = xn.position_i.norm();
    if (!(r_norm > 0.0)) throw std::domain_error("position norm must be positive");
    const Vec3 gravity_i = -p.gravitational_parameter_m3ps2 * xn.position_i / std::pow(r_norm, 3);
    const double altitude = r_norm - p.planet_radius_m;
    const Environment env = environment_at(p, altitude);
    const Vec3 rel_v_i = xn.velocity_i - env.wind_i_mps;
    const Vec3 rel_v_b = r_bi.transpose() * rel_v_i;
    const double speed = rel_v_b.norm();
    const double mach = speed / std::max(env.speed_of_sound_mps, 1e-9);

    Vec3 force_b = Vec3::Zero();
    Vec3 moment_b = Vec3::Zero();
    if (speed > 1e-9 && env.pressure_pa > 0.0) {
        const double qbar = 0.5 * env.density_kgpm3 * speed * speed;
        const double cd = std::max(0.0, interpolate(p.aero.drag_vs_mach, mach));
        const double alpha = std::atan2(std::sqrt(rel_v_b.y() * rel_v_b.y() + rel_v_b.z() * rel_v_b.z()),
                                        std::max(std::abs(rel_v_b.x()), 1e-9));
        const double cl = interpolate(p.aero.lift_vs_alpha, alpha);
        const Vec3 drag_b = -qbar * p.reference_area_m2 * cd * rel_v_b.normalized();
        // Pick a deterministic body reference axis for the lift plane. The
        // fallback avoids a 0/0 when velocity is parallel to the primary axis.
        Vec3 lift_plane_normal = rel_v_b.cross(Vec3::UnitY());
        if (lift_plane_normal.squaredNorm() < 1e-18)
            lift_plane_normal = rel_v_b.cross(Vec3::UnitZ());
        const Vec3 lift_dir = lift_plane_normal.cross(rel_v_b).normalized();
        force_b += drag_b + qbar * p.reference_area_m2 * cl * lift_dir;
    }

    double mdot = 0.0;
    if (!p.motor.empty() && xn.propellant_mass > 0.0) {
        const double thrust_vac = motor_value(p.motor, u.time_s, true);
        mdot = std::max(0.0, motor_value(p.motor, u.time_s, false));
        const double thrust = thrust_vac - p.nozzle_exit_area_m2 * env.pressure_pa;
        const Eigen::AngleAxisd tvc_yaw(u.tvc_yaw_rad + p.thrust_misalignment_yaw_rad, Vec3::UnitZ());
        const Eigen::AngleAxisd tvc_pitch(u.tvc_pitch_rad + p.thrust_misalignment_pitch_rad, Vec3::UnitY());
        const Vec3 thrust_force_b = tvc_yaw * tvc_pitch * Vec3(std::max(0.0, thrust), 0.0, 0.0);
        force_b += thrust_force_b;
        // TVC produces a rigid-body moment when the nozzle line is offset from CG.
        moment_b += p.thrust_application_point_b_m.cross(thrust_force_b);
    } else {
        mdot = 0.0;
    }
    if (xn.propellant_mass <= 0.0) mdot = 0.0;

    const Mat3 I = inertia_body(xn, p);
    StateDerivative d;
    d.position_i_dot = xn.velocity_i;
    d.velocity_i_dot = gravity_i + r_bi * force_b / mass(xn, p);
    const Eigen::Quaterniond omega_q(0.0, xn.angular_rate_b.x(), xn.angular_rate_b.y(), xn.angular_rate_b.z());
    d.attitude_bi_dot.coeffs() = (xn.attitude_bi * omega_q).coeffs() * 0.5;
    d.angular_rate_b_dot = I.inverse() * (moment_b - xn.angular_rate_b.cross(I * xn.angular_rate_b));
    d.propellant_mass_dot = -mdot;
    return d;
}

State rk4_step(const State& x, const VehicleParameters& p, const Inputs& u, double dt_s) {
    if (!(dt_s > 0.0)) throw std::invalid_argument("RK4 timestep must be positive");
    const StateDerivative k1 = derivative(x, p, u);
    const StateDerivative k2 = derivative(add_scaled(x, k1, 0.5 * dt_s), p, u);
    const StateDerivative k3 = derivative(add_scaled(x, k2, 0.5 * dt_s), p, u);
    const StateDerivative k4 = derivative(add_scaled(x, k3, dt_s), p, u);
    State y = add_scaled(x, k1, dt_s / 6.0);
    y = add_scaled(y, k2, dt_s / 3.0);
    y = add_scaled(y, k3, dt_s / 3.0);
    y = add_scaled(y, k4, dt_s / 6.0);
    normalize(y);
    y.propellant_mass = std::clamp(y.propellant_mass, 0.0, p.initial_propellant_mass_kg);
    return y;
}

double specific_orbital_energy(const State& x, const VehicleParameters& p) {
    return 0.5 * x.velocity_i.squaredNorm() - p.gravitational_parameter_m3ps2 / x.position_i.norm();
}

}  // namespace lv::dynamics
