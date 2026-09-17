#include "sensors/sensors.hpp"

namespace lv::sensors {
void MeasurementQueue::push(const Measurement& measurement) {
    double arrival=0.0;
    if (const auto* gps=std::get_if<GpsMeasurement>(&measurement)) arrival=gps->arrival_time_s;
    else if (const auto* imu=std::get_if<ImuMeasurement>(&measurement)) arrival=imu->time_s;
    else arrival=std::get<BaroMeasurement>(measurement).time_s;
    queue_.push({arrival,measurement});
}
std::vector<Measurement> MeasurementQueue::release(double current_time_s) {
    std::vector<Measurement> released;
    while(!queue_.empty() && queue_.top().arrival_time_s <= current_time_s + 1e-12) {
        released.push_back(queue_.top().measurement); queue_.pop();
    }
    return released;
}
} // namespace lv::sensors
