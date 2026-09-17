cmake_minimum_required(VERSION 3.20)
project(launch_vehicle_simulator LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(GTest REQUIRED)

add_library(lv_dynamics
    src/dynamics/dynamics.cpp
    src/sensors/sensors.cpp
    src/sensors/measurement_queue.cpp
    src/estimator/ukf.cpp
    src/control/control.cpp
    src/monte_carlo/monte_carlo.cpp
)
target_include_directories(lv_dynamics PUBLIC include)
target_link_libraries(lv_dynamics PUBLIC Eigen3::Eigen)
target_compile_options(lv_dynamics PRIVATE -Wall -Wextra -Wpedantic)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(lv_dynamics PUBLIC OpenMP::OpenMP_CXX)
endif()

enable_testing()
add_executable(dynamics_tests tests/dynamics_test.cpp tests/sensors_test.cpp tests/measurement_queue_test.cpp tests/control_test.cpp tests/estimator_test.cpp)
target_link_libraries(dynamics_tests PRIVATE lv_dynamics GTest::gtest_main)
include(GoogleTest)
gtest_discover_tests(dynamics_tests)

add_executable(run_monte_carlo src/monte_carlo/main.cpp)
target_link_libraries(run_monte_carlo PRIVATE lv_dynamics)
