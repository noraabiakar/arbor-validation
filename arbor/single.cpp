/*
 * A miniapp that demonstrates how to make a ring model
 *
 */

#include <fstream>
#include <iomanip>
#include <iostream>

#include <arbor/assert_macro.hpp>
#include <arbor/common_types.hpp>
#include <arbor/context.hpp>
#include <arbor/load_balance.hpp>
#include <arbor/cable_cell.hpp>
#include <arbor/profile/meter_manager.hpp>
#include <arbor/profile/profiler.hpp>
#include <arbor/simple_sampler.hpp>
#include <arbor/simulation.hpp>
#include <arbor/recipe.hpp>
#include <arbor/version.hpp>

#include <arborenv/concurrency.hpp>
#include <arborenv/gpu_env.hpp>

#ifdef ARB_MPI_ENABLED
#include <mpi.h>
#include <arborenv/with_mpi.hpp>
#endif

#include "parameters.hpp"

using arb::cell_gid_type;
using arb::cell_lid_type;
using arb::cell_size_type;
using arb::cell_member_type;
using arb::cell_kind;
using arb::time_type;
using arb::cell_probe_address;

// Writes voltage trace as a json file.
void write_trace_json(const arb::trace_data<double>& trace);

// Generate a cell.
arb::cable_cell single_cell(const single_params& params);

class soma_recipe: public arb::recipe {
public:
    soma_recipe(single_params params):
      num_cells_(1), params_(params), event_weight_(params.weight), catalogue_(arb::global_default_catalogue()) {
        cell_gprop_.catalogue = &catalogue_;
        cell_gprop_.default_parameters = arb::neuron_parameter_defaults;
        cell_gprop_.default_parameters.temperature_K = params_.temp + 273.15;
        cell_gprop_.default_parameters.init_membrane_potential = params_.v_init;

        if (params.mech == "borgka" || params.mech == "cagk") {
            cell_gprop_.default_parameters.ion_data["k"].init_reversal_potential = -77.0;
        }
    }

    cell_size_type num_cells() const override {
        return num_cells_;
    }

    arb::util::unique_any get_cell_description(cell_gid_type gid) const override {
        return single_cell(params_);
    }

    cell_kind get_cell_kind(cell_gid_type gid) const override {
        return cell_kind::cable;
    }

    // Each cell has one spike detector (at the soma).
    cell_size_type num_sources(cell_gid_type gid) const override {
        return 0;
    }

    cell_size_type num_targets(cell_gid_type gid) const override {
        return 1;
    }

    // Return one event generator on gid 0. This generates a single event that will
    // kick start the spiking.
    std::vector<arb::event_generator> event_generators(cell_gid_type gid) const override {
        std::vector<arb::event_generator> gens;
        arb::pse_vector svec;

        std::vector<double> spikes = {
                25.269724183039855, 29.37076391451496,
                58.472477010286546, 93.80268485203328,
                112.71090127018375, 142.6472406502223
        };

        for (auto s: spikes) {
            svec.push_back({{gid, 0}, s, event_weight_});
        }
        gens.push_back(arb::explicit_generator(svec));
        return gens;
    }

    // There is one probe (for measuring voltage at the soma) on the cell.
    cell_size_type num_probes(cell_gid_type gid)  const override {
        return 1;
    }

    arb::probe_info get_probe(cell_member_type id) const override {
        // Get the appropriate kind for measuring voltage.
        cell_probe_address::probe_kind kind = cell_probe_address::membrane_voltage;
        // Measure at the soma.
        arb::segment_location loc(0, 0.5);

        return arb::probe_info{id, kind, cell_probe_address{loc, kind}};
    }

    arb::util::any get_global_properties(cell_kind k) const override {
        return cell_gprop_;
    }

    void add_ion(const std::string& ion_name, int charge, double init_iconc, double init_econc, double init_revpot) {
        cell_gprop_.add_ion(ion_name, charge, init_iconc, init_econc, init_revpot);
    }


private:
    cell_size_type num_cells_;
    single_params params_;
    float event_weight_;

    arb::cable_cell_global_properties cell_gprop_;
    arb::mechanism_catalogue catalogue_;
};


int main(int argc, char** argv) {
    try {
        bool root = true;

        arb::proc_allocation resources;
        if (auto nt = arbenv::get_env_num_threads()) {
            resources.num_threads = nt;
        }
        else {
            resources.num_threads = arbenv::thread_concurrency();
        }

#ifdef ARB_MPI_ENABLED
        arbenv::with_mpi guard(argc, argv, false);
        resources.gpu_id = arbenv::find_private_gpu(MPI_COMM_WORLD);
        auto context = arb::make_context(resources, MPI_COMM_WORLD);
        root = arb::rank(context) == 0;
#else
        resources.gpu_id = arbenv::default_gpu();
        auto context = arb::make_context(resources);
#endif

#ifdef ARB_PROFILE_ENABLED
        arb::profile::profiler_initialize(context);
#endif

        // Print a banner with information about hardware configuration
        std::cout << "gpu:      " << (has_gpu(context)? "yes": "no") << "\n";
        std::cout << "threads:  " << num_threads(context) << "\n";
        std::cout << "mpi:      " << (has_mpi(context)? "yes": "no") << "\n";
        std::cout << "ranks:    " << num_ranks(context) << "\n" << std::endl;

        arb::profile::meter_manager meters;
        meters.start(context);

        // Create an instance of our recipe.
        auto params = read_params(argc, argv);
        soma_recipe recipe(params);

        if(params.mech == "cagk" || params.mech == "ccanl") {
            recipe.add_ion("nca", 2, 1, 2.0/3, 0);
            recipe.add_ion("lca", 2, 1, 2.0/3, 0);
            recipe.add_ion("tca", 2, 1, 2.0/3, 0);
        }
        else if(params.mech == "cat") {
            recipe.add_ion("tca", 2, 0, 2.0/3, 0);
        }
        else if(params.mech == "lca") {
            recipe.add_ion("lca", 2, 0, 2.0/3, 0);
        }
        else if(params.mech == "nca") {
            recipe.add_ion("nca", 2, 0, 2.0/3, 0);
        }
        else if(params.mech == "gskch") {
            recipe.add_ion("nca", 2, 0, 2.0/3, 0);
            recipe.add_ion("lca", 2, 0, 2.0/3, 0);
            recipe.add_ion("tca", 2, 0, 2.0/3, 0);
            recipe.add_ion("sk", 1, 0, 2.0/3, 0);
        }
        if(params.mech == "ichan2") {
            recipe.add_ion("nat", 1, 0, 2.0/3, 0);
            recipe.add_ion("kf", 1, 0, 2.0/3, 0);
            recipe.add_ion("ks", 1, 0, 2.0/3, 0);
        }

        auto decomp = arb::partition_load_balance(recipe, context);

        // Construct the model.
        arb::simulation sim(recipe, decomp, context);

        // Set up the probe that will measure voltage in the cell.

        // The id of the only probe on the cell: the cell_member type points to (cell 0, probe 0)
        auto probe_id = cell_member_type{0, 0};

        // The schedule for sampling is 10 samples every 1 ms.
        auto sched = arb::regular_schedule(0.001);

        // This is where the voltage samples will be stored as (time, value) pairs
        arb::trace_data<double> voltage;
        // Now attach the sampler at probe_id, with sampling schedule sched, writing to voltage
        sim.add_sampler(arb::one_probe(probe_id), sched, arb::make_simple_sampler(voltage));

        // Set up recording of spikes to a vector on the root process.
        std::vector<arb::spike> recorded_spikes;
        if (root) {
            sim.set_global_spike_callback(
                [&recorded_spikes](const std::vector<arb::spike>& spikes) {
                    recorded_spikes.insert(recorded_spikes.end(), spikes.begin(), spikes.end());
                });
        }

        meters.checkpoint("model-init", context);

        std::cout << "running simulation" << std::endl;
        // Run the simulation for 100 ms, with time steps of 0.025 ms.
        sim.run(200, params.dt);

        meters.checkpoint("model-run", context);

        auto ns = sim.num_spikes();

        // Write spikes to file
        if (root) {
            std::cout << "\n" << ns << " spikes generated\n.";
            std::ofstream fid("spikes.gdf");
            if (!fid.good()) {
                std::cerr << "Warning: unable to open file spikes.gdf for spike output\n";
            }
            else {
                char linebuf[45];
                for (auto spike: recorded_spikes) {
                    auto n = std::snprintf(
                        linebuf, sizeof(linebuf), "%u %.4f\n",
                        unsigned{spike.source.gid}, float(spike.time));
                    fid.write(linebuf, n);
                }
            }
        }

        // Write the samples to a json file.
        if (root) write_trace_json(voltage);

        auto report = arb::profile::make_meter_report(meters, context);
        std::cout << report;
    }
    catch (std::exception& e) {
        std::cerr << "exception caught in ring miniapp: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

void write_trace_json(const arb::trace_data<double>& trace) {
    std::string path = "./voltages.json";

    nlohmann::json json;
    json["name"] = "ring demo";
    json["units"] = "mV";
    json["cell"] = "0.0";
    json["probe"] = "0";

    auto& jt = json["data"]["time"];
    auto& jy = json["data"]["voltage"];

    for (const auto& sample: trace) {
        jt.push_back(sample.t);
        jy.push_back(sample.v);
    }

    std::ofstream file(path);
    file << std::setw(1) << json << "\n";
}

arb::cable_cell single_cell(const single_params& params) {
    arb::cable_cell cell;

    // Add soma.
    auto soma = cell.add_soma(11.65968/2.0);

    auto dend = cell.add_cable(0, arb::section_kind::dendrite, 30.0/2.0, 30.0/2.0, 200);
    dend->set_compartments(2000);

    if (params.soma_mech) {
        auto mech = arb::mechanism_desc(params.mech);
        if (params.mech == "borgka") {
            mech.set("gkabar", params.gkabar);
        }
        soma->add_mechanism(mech);

        if (params.mech == "ccanl") {
            cell.default_parameters.reversal_potential_method["nca"] = "ccanlrev";
            cell.default_parameters.reversal_potential_method["lca"] = "ccanlrev";
            cell.default_parameters.reversal_potential_method["tca"] = "ccanlrev";
        }

    } else {
        auto pas = arb::mechanism_desc("pas");
        pas.set("g", params.pas_g);
        pas.set("e", params.pas_e);

        soma->add_mechanism(pas);
    }

    if (params.dend_mech) {
        auto mech = arb::mechanism_desc(params.mech);
        if (params.mech == "borgka") {
            mech.set("gkabar", params.gkabar);
        }
        if (params.mech == "ccanl") {
            cell.default_parameters.reversal_potential_method["nca"] = "ccanlrev";
            cell.default_parameters.reversal_potential_method["lca"] = "ccanlrev";
            cell.default_parameters.reversal_potential_method["tca"] = "ccanlrev";
        }

        dend->add_mechanism(mech);
    } else {
        auto pas = arb::mechanism_desc("pas");
        pas.set("g", params.pas_g);
        pas.set("e", params.pas_e);

        dend->add_mechanism(pas);
    }

    auto exp2syn = arb::mechanism_desc("exp2syn");
    exp2syn.set("tau1", params.tau1_syn);
    exp2syn.set("tau2", params.tau2_syn);
    exp2syn.set("e", params.e_syn);

    cell.add_synapse({params.syn_seg, params.syn_loc}, exp2syn);
    std::cout << params.syn_seg << " " << params.syn_loc << std::endl;

    return cell;
}

