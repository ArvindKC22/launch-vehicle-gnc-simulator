#pragma once
#include "control/control.hpp"
#include <string>

namespace lv::monte_carlo {
struct RunResult { int run; double apogee_m; double impact_x_m; double impact_y_m;
    double wind_scale; double thrust_scale; double dry_mass_scale; double aero_scale; double sensor_bias_scale;
    double thrust_misalignment_deg; double cg_offset_m; int impact_valid; };
struct Config { int runs{10000}; double duration_s{300.0}; double dt_s{0.01}; unsigned threads{1}; std::string output_csv{"mc_results.csv"}; };
std::vector<RunResult> run(const dynamics::VehicleParameters& vehicle, const sensors::SensorConfig& sensor,
                           const Config& config);
void write_csv(const std::vector<RunResult>& results, const std::string& path);
}
