#include "sensors/sensors.hpp"
#include <gtest/gtest.h>
TEST(Sensors, RatesAndLatencyAreExplicit) {
    lv::dynamics::VehicleParameters p; p.dry_mass_kg=10; p.planet_radius_m=6378137;
    lv::dynamics::State x; x.position_i={6379137,0,0}; x.propellant_mass=0;
    lv::sensors::SensorSuite s({},7); auto dx=lv::dynamics::derivative(x,p,{});
    auto z=s.sample(x,dx,0.02,0.02,p); int imu_count=0; bool gps=false,baro=false;
    for(auto& m:z){imu_count+=std::holds_alternative<lv::sensors::ImuMeasurement>(m); gps|=std::holds_alternative<lv::sensors::GpsMeasurement>(m); baro|=std::holds_alternative<lv::sensors::BaroMeasurement>(m);}
    EXPECT_EQ(imu_count,5); EXPECT_TRUE(gps); EXPECT_TRUE(baro);
    for(auto& m:z) if(auto* g=std::get_if<lv::sensors::GpsMeasurement>(&m)) EXPECT_DOUBLE_EQ(g->arrival_time_s,g->measurement_time_s+.2);
}

TEST(Sensors, TurnOnBiasIsDrawnOnce) {
    lv::sensors::SensorConfig c;
    lv::sensors::SensorSuite s(c,42);
    EXPECT_GT(s.accel_bias().norm(),0.0);
    EXPECT_GT(s.gyro_bias().norm(),0.0);
}

TEST(Sensors, BarometerLagUsesBarometerSamplePeriod) {
    lv::dynamics::VehicleParameters p; p.dry_mass_kg=10; p.planet_radius_m=6378137;
    lv::dynamics::State x; x.position_i={p.planet_radius_m+1000.0,0,0};
    lv::sensors::SensorConfig c; c.baro_pressure_sigma_pa=0.0;
    lv::sensors::SensorSuite s(c,7); const auto dx=lv::dynamics::derivative(x,p,{});
    const auto measurements=s.sample(x,dx,0.0,0.002,p);
    double altitude=-1.0;
    for(const auto& m:measurements)
        if(const auto* baro=std::get_if<lv::sensors::BaroMeasurement>(&m)) altitude=baro->altitude_m;
    const double expected=(1.0/c.baro_rate_hz)/(c.baro_lag_s+1.0/c.baro_rate_hz)*1000.0;
    EXPECT_NEAR(altitude,expected,0.5);
}
