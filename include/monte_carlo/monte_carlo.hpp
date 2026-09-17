#pragma once
#include "control/control.hpp"
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace lv::monte_carlo {
struct RunResult { int run; double apogee_m; double impact_x_m; double impact_y_m;
    double wind_scale; double thrust_scale; double dry_mass_scale; double aero_scale; double sensor_bias_scale;
    double thrust_misalignment_deg; double cg_offset_m; double mean_nees; double mean_nis;
    double mean_gps_nis; double mean_baro_nis; int impact_valid; };
struct ConsistencySample {
    int run;
    double time_s;
    std::string metric;
    double value;
    int dof;
};
struct Config {
    int runs{10000};
    double duration_s{300.0};
    double dt_s{0.01};
    unsigned threads{1};
    bool closed_loop{false};
    bool navigation_only{false};
    bool truth_guidance{false};
    double consistency_interval_s{0.5};
    std::string output_csv{"mc_results.csv"};
    std::string consistency_output_csv{};
    std::uint64_t seed{20260908};
    double wind_sigma{0.20};
    double thrust_sigma{0.02};
    double dry_mass_sigma{0.01};
    double aero_sigma{0.03};
    double sensor_bias_sigma{0.20};
    double thrust_misalignment_sigma_deg{0.15};
    double cg_offset_sigma_m{0.001};
};
std::vector<RunResult> run(const dynamics::VehicleParameters& vehicle, const sensors::SensorConfig& sensor,
                           const Config& config,
                           std::vector<ConsistencySample>* consistency_samples=nullptr);
void write_csv(const std::vector<RunResult>& results, const std::string& path);
void write_consistency_csv(const std::vector<ConsistencySample>& samples, const std::string& path);
}
