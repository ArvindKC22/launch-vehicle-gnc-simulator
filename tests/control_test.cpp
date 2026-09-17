#include "control/control.hpp"
#include <gtest/gtest.h>
TEST(Control, ActuatorHonorsSlewAndSaturation) {
    lv::control::TvcrActuator a({0.1,0.5,0.08}); auto y=a.update({1,1},.01);
    EXPECT_LE(y.cwiseAbs().maxCoeff(),.08); EXPECT_LE(y.cwiseAbs().maxCoeff(),.005+1e-12);
}
