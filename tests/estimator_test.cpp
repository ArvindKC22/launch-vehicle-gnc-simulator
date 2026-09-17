#include "estimator/ukf.hpp"
#include <gtest/gtest.h>
TEST(Estimator, MRPErrorDoesNotAverageQuaternions) {
    lv::dynamics::VehicleParameters p; p.dry_mass_kg=10; p.planet_radius_m=6378137;
    lv::dynamics::State x; x.position_i={6379137,0,0}; x.propellant_mass=0;
    Eigen::Matrix<double,lv::estimator::kErrorN,lv::estimator::kErrorN> P=Eigen::Matrix<double,lv::estimator::kErrorN,lv::estimator::kErrorN>::Identity()*1e-4;
    lv::estimator::UKF f(x,P); f.propagate(.01,p,{}); EXPECT_NEAR(f.estimate().nominal.attitude_bi.norm(),1.0,1e-14);
}

TEST(Estimator, EKFHasIndependentPropagationAndUpdateInterface) {
    lv::dynamics::VehicleParameters p; p.dry_mass_kg=10; p.planet_radius_m=6378137;
    lv::dynamics::State x; x.position_i={6379137,0,0}; x.propellant_mass=0;
    Eigen::Matrix<double,lv::estimator::kErrorN,lv::estimator::kErrorN> P=Eigen::Matrix<double,lv::estimator::kErrorN,lv::estimator::kErrorN>::Identity()*1e-4;
    lv::estimator::EKF f(x,P); f.propagate(.01,p,{}); EXPECT_NEAR(f.estimate().nominal.attitude_bi.norm(),1.0,1e-14);
}

TEST(Estimator, UKFAcceptsAsynchronousImuRateMeasurement) {
    lv::dynamics::VehicleParameters p; p.dry_mass_kg=10; p.planet_radius_m=6378137;
    lv::dynamics::State x; x.position_i={6379137,0,0}; x.propellant_mass=0;
    Eigen::Matrix<double,lv::estimator::kErrorN,lv::estimator::kErrorN> P=Eigen::Matrix<double,lv::estimator::kErrorN,lv::estimator::kErrorN>::Identity()*1e-4;
    lv::estimator::UKF f(x,P); lv::sensors::ImuMeasurement imu{0.005,Eigen::Vector3d::Zero(),Eigen::Vector3d(0.1,0.2,0.3)}; f.update(imu,p);
    EXPECT_EQ(f.last_nis_dof(),3);
    EXPECT_TRUE(std::isfinite(f.last_nis()));
    EXPECT_GT(f.estimate().gyro_bias.norm(),0.0);
}

TEST(Estimator, ImuPropagationUsesSpecificForceForNavigation) {
    lv::dynamics::VehicleParameters p;
    p.dry_mass_kg=10.0; p.planet_radius_m=6378137.0;
    p.gravitational_parameter_m3ps2=0.0;
    lv::dynamics::State x; x.position_i={6378137.0,0.0,0.0};
    Eigen::Matrix<double,lv::estimator::kErrorN,lv::estimator::kErrorN> P=
        Eigen::Matrix<double,lv::estimator::kErrorN,lv::estimator::kErrorN>::Identity()*1e-12;
    lv::estimator::UKF f(x,P);
    f.update(lv::sensors::ImuMeasurement{0.0,Eigen::Vector3d(1.0,0.0,0.0),
                                         Eigen::Vector3d::Zero()},p);
    f.propagate(0.1,p,{});
    EXPECT_NEAR(f.estimate().nominal.velocity_i.x(),0.1,1e-6);
    EXPECT_NEAR(f.estimate().nominal.position_i.x(),6378137.005,1e-6);
}

TEST(Estimator, GpsLatencyIsCompensatedToFilterTime) {
    lv::dynamics::VehicleParameters p; p.dry_mass_kg=10; p.planet_radius_m=6378137; p.gravitational_parameter_m3ps2=0.0;
    lv::dynamics::State x; x.position_i={6379137,0,0}; x.velocity_i={0,10,0};
    Eigen::Matrix<double,lv::estimator::kErrorN,lv::estimator::kErrorN> P=
        Eigen::Matrix<double,lv::estimator::kErrorN,lv::estimator::kErrorN>::Identity();
    lv::estimator::UKF f(x,P);
    f.propagate(0.2,p,{});
    const auto predicted=f.estimate().nominal;
    lv::sensors::GpsMeasurement gps{0.0,0.2,x.position_i,x.velocity_i};
    f.update(gps,p);
    EXPECT_LT((f.estimate().nominal.position_i-predicted.position_i).norm(),1e-8);
    EXPECT_EQ(f.last_nis_dof(),6);
}

TEST(Estimator, MRPShadowSetAvoidsThreeSixtyDegreeSingularity) {
    Eigen::Vector3d sigma(2.0,0.0,0.0); auto shadow=lv::estimator::mrp_shadow_set(sigma);
    EXPECT_NEAR(shadow.norm(),0.5,1e-15); EXPECT_LT(shadow.norm(),1.0);
}

TEST(Estimator, LinearUkfMatchesAnalyticKalmanFilter) {
    Eigen::VectorXd x(2); x<<10.0,2.0; Eigen::MatrixXd P=Eigen::Matrix2d::Identity();
    Eigen::MatrixXd F(2,2); F<<1,1,0,1; Eigen::MatrixXd Q=Eigen::Matrix2d::Identity()*0.01;
    Eigen::MatrixXd H(1,2); H<<1,0; Eigen::MatrixXd R(1,1); R<<0.25; Eigen::VectorXd z(1); z<<12.4;
    auto ukf=lv::estimator::linear_ukf_step(x,P,F,Q,H,R,z);
    auto xp=F*x; auto Pp=F*P*F.transpose()+Q; auto S=H*Pp*H.transpose()+R; auto K=Pp*H.transpose()*S.inverse(); auto kf=xp+K*(z-H*xp); auto Pkf=Pp-K*S*K.transpose();
    // The two paths perform different floating-point operation orderings;
    // 1e-10 is tight relative agreement for this double-precision benchmark.
    EXPECT_LT((ukf.state-kf).norm(),1e-10); EXPECT_LT((ukf.covariance-Pkf).norm(),1e-10);
}
