#include <iostream>

#include <array>
#include <cmath>
#include <fstream>
#include <random>

#include <arbor/cable_cell.hpp>
#include <common/json_params.hpp>

std::vector<double> read_spike_times();

struct single_params {
    double temp, v_init;
    double tau1_syn, tau2_syn, e_syn;
    double hh_gnabar, hh_gkbar, hh_gl, hh_ena, hh_ek;
    double pas_e, pas_g;
    unsigned syn_seg;
    double syn_loc;
    double dt, weight;
    bool soma_mech, dend_mech;
    std::string mech;
};

single_params read_params(int argc, char** argv) {
    single_params p;

    using sup::param_from_json;

    if (argc<2) {
        throw std::runtime_error("No input parameter file provided.");
    }
    if (argc>2) {
        throw std::runtime_error("More than command line one option not permitted.");
    }

    std::string fname = argv[1];
    std::cout << "Loading parameters from file: " << fname << "\n";
    std::ifstream f(fname);

    if (!f.good()) {
        throw std::runtime_error("Unable to open input parameter file: "+fname);
    }

    nlohmann::json json;
    json << f;

    param_from_json(p.temp, "temp", json);
    param_from_json(p.v_init, "vinit", json);
    param_from_json(p.dt, "dt_arbor", json);
    param_from_json(p.tau1_syn, "tau1_syn", json);
    param_from_json(p.tau2_syn, "tau2_syn", json);
    param_from_json(p.e_syn, "e_syn", json);
    param_from_json(p.hh_gnabar, "hh_gnabar", json);
    param_from_json(p.hh_gkbar, "hh_gkbar", json);
    param_from_json(p.hh_gl, "hh_gl", json);
    param_from_json(p.hh_ena, "hh_ena", json);
    param_from_json(p.hh_ek, "hh_ek", json);
    param_from_json(p.pas_e, "pas_e", json);
    param_from_json(p.pas_g, "pas_g", json);
    param_from_json(p.syn_seg, "syn_seg", json);
    param_from_json(p.syn_loc, "syn_loc", json);
    param_from_json(p.weight, "weight", json);
    param_from_json(p.soma_mech, "soma_mech", json);
    param_from_json(p.dend_mech, "dend_mech", json);
    param_from_json(p.mech, "mech", json);

    for (auto it=json.begin(); it!=json.end(); ++it) {
        std::cout << "  Warning: unused input parameter: \"" << it.key() << "\"\n";
    }
    std::cout << "\n";

    return p;

}