#include "dynamics/dynamics.hpp"

#include <gtest/gtest.h>

namespace {
using namespace lv::dynamics;

VehicleParameters ballistic_vehicle() {
    VehicleParameters p;
    p.dry_mass_kg = 100.0;
    p.initial_propellant_mass_kg = 0.0;
    p.reference_area_m2 = 0.01;
    p.planet_radius_m = 6378137.0;
    p.aero.drag_vs_mach = {{0.0, 0.0}, {20.0, 0.0}};
    return p;
}

TEST(Dynamics, DragFreeUnpoweredArcConservesSpecificOrbitalEnergy) {
    const auto p = ballistic_vehicle();
    State x;
    x.position_i = Vec3(p.planet_radius_m + 120000.0, 0.0, 0.0);
    x.velocity_i = Vec3(0.0, 7700.0, 0.0);
    x.propellant_mass = 0.0;
    const double e0 = specific_orbital_energy(x, p);

    // dt=0.1 s and RK4's O(dt^4) global truncation error define the numerical
    // tolerance here. The bound is 1e-9 relative energy over a 600 s arc.
    for (int i = 0; i < 6000; ++i) x = rk4_step(x, p, Inputs{static_cast<double>(i) * 0.1, 0.0, 0.0}, 0.1);
    const double e1 = specific_orbital_energy(x, p);
    EXPECT_LT(std::abs((e1 - e0) / e0), 1e-9);
}

TEST(Dynamics, Rk4NormalizesQuaternion) {
    auto p = ballistic_vehicle();
    State x;
    x.position_i = Vec3(p.planet_radius_m + 1000.0, 0.0, 0.0);
    x.velocity_i = Vec3(0.0, 10.0, 0.0);
    x.angular_rate_b = Vec3(0.1, -0.2, 0.3);
    x = rk4_step(x, p, Inputs{}, 0.5);
    EXPECT_NEAR(x.attitude_bi.norm(), 1.0, 1e-15);
}

}  // namespace
