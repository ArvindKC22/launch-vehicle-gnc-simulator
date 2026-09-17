#include "monte_carlo/monte_carlo.hpp"
#include <fstream>
#include <regex>
#include <sstream>
#include <iostream>
#include <thread>
#ifdef _OPENMP
#include <omp.h>
#endif
static std::vector<double> numbers(const std::string& s) { std::vector<double> v; std::regex r(R"([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)"); for(std::sregex_iterator i(s.begin(),s.end(),r),e;i!=e;++i)v.push_back(std::stod(i->str())); return v; }
static std::string line_for(const std::string& path,const std::string& key) { std::ifstream f(path); std::string l; while(std::getline(f,l)) if(l.find(key)!=std::string::npos)return l; return {}; }
static lv::dynamics::VehicleParameters load_vehicle(const std::string& path) {
    using namespace lv::dynamics; VehicleParameters p; auto scalar=[&](const char* k,double d){auto n=numbers(line_for(path,k));return n.empty()?d:n.back();};
    p.dry_mass_kg=scalar("dry_mass_kg",100); p.initial_propellant_mass_kg=scalar("initial_propellant_mass_kg",50); p.reference_area_m2=scalar("reference_area_m2",.01); p.nozzle_exit_area_m2=scalar("nozzle_exit_area_m2",.0002);
    auto a=numbers(line_for(path,"drag_vs_mach")); for(size_t i=0;i+1<a.size();i+=2)p.aero.drag_vs_mach.push_back({a[i],a[i+1]});
    a=numbers(line_for(path,"lift_vs_alpha")); for(size_t i=0;i+1<a.size();i+=2)p.aero.lift_vs_alpha.push_back({a[i],a[i+1]});
    a=numbers(line_for(path,"motor:")); for(size_t i=0;i+2<a.size();i+=3)p.motor.push_back({a[i],a[i+1],a[i+2]});
    a=numbers(line_for(path,"wind:")); for(size_t i=0;i+2<a.size();i+=3)p.wind.push_back({a[i],lv::dynamics::Vec3(0,a[i+1],a[i+2])});
    auto dry=numbers(line_for(path,"inertia_dry_b:")); if(dry.size()>=3) p.inertia_dry_b=lv::dynamics::Vec3(dry[0],dry[1],dry[2]).asDiagonal();
    auto full=numbers(line_for(path,"inertia_full_b:")); if(full.size()>=3) p.inertia_full_b=lv::dynamics::Vec3(full[0],full[1],full[2]).asDiagonal();
    auto arm=numbers(line_for(path,"thrust_application_point_b_m:")); if(arm.size()>=3) p.thrust_application_point_b_m=lv::dynamics::Vec3(arm[0],arm[1],arm[2]);
    return p;
}
int main(int argc,char** argv) {
    const std::string root=argc>2?argv[2]:"."; auto p=load_vehicle(root+"/config/vehicle.yaml");
    lv::monte_carlo::Config c; if(argc>1)c.output_csv=argv[1]; c.threads=std::max(1u,std::thread::hardware_concurrency());
    const std::string dispersion_path=root+"/config/dispersions.yaml";
    auto dscalar=[&](const char* key,double fallback){ auto values=numbers(line_for(dispersion_path,key)); return values.empty()?fallback:values.back(); };
    c.runs=static_cast<int>(dscalar("runs:",c.runs));
    c.seed=static_cast<std::uint64_t>(dscalar("seed:",c.seed));
    c.wind_sigma=dscalar("wind_magnitude_scale:",c.wind_sigma);
    c.thrust_sigma=dscalar("thrust_scale:",c.thrust_sigma);
    c.thrust_misalignment_sigma_deg=dscalar("thrust_misalignment_deg:",c.thrust_misalignment_sigma_deg);
    c.dry_mass_sigma=dscalar("dry_mass_scale:",c.dry_mass_sigma);
    c.cg_offset_sigma_m=dscalar("cg_location_m:",c.cg_offset_sigma_m);
    c.aero_sigma=dscalar("aero_scale:",c.aero_sigma);
    c.sensor_bias_sigma=dscalar("sensor_bias_scale:",c.sensor_bias_sigma);
    c.duration_s=dscalar("duration_s:",c.duration_s);
    c.dt_s=dscalar("dt_s:",c.dt_s);
    c.threads=static_cast<unsigned>(std::max(1.0,dscalar("threads:",c.threads)));
    c.closed_loop=dscalar("closed_loop:",c.closed_loop?1.0:0.0)!=0.0;
    c.consistency_interval_s=dscalar("consistency_interval_s:",c.consistency_interval_s);
    // Optional: output root runs duration_s dt_s threads. This makes smoke
    // tests practical without editing the declarative dispersion file.
    if(argc>3)c.runs=std::stoi(argv[3]); if(argc>4)c.duration_s=std::stod(argv[4]); if(argc>5)c.dt_s=std::stod(argv[5]); if(argc>6)c.threads=static_cast<unsigned>(std::stoul(argv[6])); if(argc>7)c.closed_loop=(std::stoi(argv[7])!=0); if(argc>8)c.consistency_output_csv=argv[8]; if(argc>9)c.consistency_interval_s=std::stod(argv[9]); if(argc>10)c.navigation_only=(std::stoi(argv[10])!=0); if(argc>11)c.truth_guidance=(std::stoi(argv[11])!=0);
#ifdef _OPENMP
    omp_set_num_threads(static_cast<int>(c.threads));
#endif
    std::vector<lv::monte_carlo::ConsistencySample> consistency;
    auto r=lv::monte_carlo::run(p,{},c,
        c.consistency_output_csv.empty()?nullptr:&consistency);
    lv::monte_carlo::write_csv(r,c.output_csv);
    std::cout<<"Wrote "<<r.size()<<" runs to "<<c.output_csv<<"\n";
    if(!c.consistency_output_csv.empty()) {
        lv::monte_carlo::write_consistency_csv(consistency,c.consistency_output_csv);
        std::cout<<"Wrote "<<consistency.size()<<" consistency samples to "
                 <<c.consistency_output_csv<<"\n";
    }
}
