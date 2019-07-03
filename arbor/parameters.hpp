#include <iostream>

#include <array>
#include <cmath>
#include <fstream>
#include <random>

#include <arbor/cable_cell.hpp>
#include <common/json_params.hpp>

std::vector<double> read_spike_times();

struct single_params {
    std::vector<double> spikes;
    double temp, v_init;
    double dt, weight;
    double tau1_syn, tau2_syn, e_syn;

    double syn_loc;

    double pas_e, pas_g;

    double gnatbar_ichan2, gkfbar_ichan2, gksbar_ichan2, gl_ichan2;
    double gkabar_borgka;
    double gncabar_nca;
    double glcabar_lca;
    double gcatbar_cat;
    double gskbar_gskch;
    double gkbar_cagk;
    double catau_ccanl, caiinf_ccanl;

    double ra_mult, cm_mult, ra, cm;

    double enat, ekf, eks, ek, elca, etca, esk, el_ichan2, cao;
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

    param_from_json(p.spikes, "spikes", json);
    param_from_json(p.temp, "temp", json);
    param_from_json(p.v_init, "vinit", json);
    param_from_json(p.dt, "dt_arbor", json);
    param_from_json(p.weight, "weight", json);
    param_from_json(p.tau1_syn, "tau1_syn", json);
    param_from_json(p.tau2_syn, "tau2_syn", json);
    param_from_json(p.e_syn, "e_syn", json);

    param_from_json(p.syn_loc, "syn_loc", json);
    param_from_json(p.pas_e, "pas_e", json);
    param_from_json(p.pas_g, "pas_g", json);

    param_from_json(p.gnatbar_ichan2, "gnatbar_ichan2", json);
    param_from_json(p.gkfbar_ichan2, "gkfbar_ichan2", json);
    param_from_json(p.gksbar_ichan2, "gksbar_ichan2", json);
    param_from_json(p.gl_ichan2, "gl_ichan2", json);

    param_from_json(p.gkabar_borgka, "gkabar_borgka", json);

    param_from_json(p.gncabar_nca, "gncabar_nca", json);

    param_from_json(p.glcabar_lca, "glcabar_lca", json);

    param_from_json(p.gcatbar_cat, "gcatbar_cat", json);

    param_from_json(p.gskbar_gskch, "gskbar_gskch", json);

    param_from_json(p.gkbar_cagk, "gkbar_cagk", json);

    param_from_json(p.catau_ccanl, "catau_ccanl", json);
    param_from_json(p.caiinf_ccanl, "caiinf_ccanl", json);

    param_from_json(p.ra, "ra", json);
    param_from_json(p.cm, "cm", json);
    param_from_json(p.ra_mult, "ra_mult", json);
    param_from_json(p.cm_mult, "cm_mult", json);

    param_from_json(p.enat, "enat", json);
    param_from_json(p.ekf, "ekf", json);
    param_from_json(p.eks, "eks", json);
    param_from_json(p.ek, "ek", json);
    param_from_json(p.elca, "elca", json);
    param_from_json(p.etca, "etca", json);
    param_from_json(p.esk, "esk", json);
    param_from_json(p.el_ichan2, "el_ichan2", json);
    param_from_json(p.cao, "cao", json);

    for (auto it=json.begin(); it!=json.end(); ++it) {
        std::cout << "  Warning: unused input parameter: \"" << it.key() << "\"\n";
    }
    std::cout << "\n";

    return p;

}