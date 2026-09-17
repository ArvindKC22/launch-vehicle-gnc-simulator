#include "sensors/sensors.hpp"
#include <gtest/gtest.h>
TEST(Sensors, MeasurementQueueReleasesGpsAtArrivalTime) {
    lv::sensors::MeasurementQueue q;
    q.push(lv::sensors::GpsMeasurement{1.0,1.2,Eigen::Vector3d::Zero(),Eigen::Vector3d::Zero()});
    EXPECT_EQ(q.release(1.19).size(),0u);
    auto ready=q.release(1.20); EXPECT_EQ(ready.size(),1u);
    EXPECT_TRUE(std::holds_alternative<lv::sensors::GpsMeasurement>(ready.front()));
}
