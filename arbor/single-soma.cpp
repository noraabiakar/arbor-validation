/*
 * A miniapp that demonstrates how to make a ring model
 *
 */

#include <fstream>
#include <iomanip>
#include <iostream>

#include <json.hpp>

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
arb::cable_cell soma_cell();

class soma_recipe: public arb::recipe {
public:
    soma_recipe(): num_cells_(1) {}

    cell_size_type num_cells() const override {
        return num_cells_;
    }

    arb::util::unique_any get_cell_description(cell_gid_type gid) const override {
        return soma_cell();
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
                25.269724183039855, 29.37076391451496, 58.472477010286546, 93.80268485203328,
                112.71090127018375, 142.6472406502223, 293.3318516075217, 456.96763867081177,
                559.0808556112847, 1091.4908947982385, 1265.7627055896965, 1286.3817213526308,
                1726.173576102434, 1744.7268859340948, 2072.358562894649, 2429.33700123744,
                2438.882187257708, 2444.85687651729, 2500.3783411639783, 2523.5435646207175,
                2633.0843793058734, 2663.3213690478333, 3081.2091382891876, 3104.234316785872,
                3209.158889778191, 3311.7494479555003, 3628.0607334944084, 3892.4078916268645,
                3905.382779351989, 3972.478937325283, 4039.5190445966164, 4275.579471624872,
                4761.533462201488, 4875.327265653268, 4946.068519184462, 5186.38947000671,
                5250.193973512949, 5405.921064743746, 6075.548637890089, 6106.605889233615,
                6392.503123563068, 6484.87209757147, 6622.667183819736, 7132.979244248274,
                7214.140067854784, 7632.383314077037, 7662.664989661292, 7663.029726732657,
                8205.521919274834, 8514.66346930178, 8998.325190823954, 9387.223218469393,
                9453.933657798472, 9544.328220467469, 9858.711284584584, 9955.045230553718,
                9956.054906300105
        };

        for (auto s: spikes) {
            svec.push_back({{0, 0}, s, event_weight_});
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
        arb::cable_cell_global_properties a;
        a.temperature_K = 308.15;
        a.init_membrane_potential_mV = -70;
        return a;
    }

private:
    cell_size_type num_cells_;
    float event_weight_ = 1.17;
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
        soma_recipe recipe;

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
        sim.run(10000, 0.0025);

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

arb::cable_cell soma_cell() {
    arb::cable_cell cell;

    // Add soma.
    auto soma = cell.add_soma(11.65968/2.0); // For area of 500 μm².

    auto hh = arb::mechanism_desc("hh");
    hh.set("ena", 50);
    hh.set("ek", -77);
    soma->add_mechanism(hh);

    auto expsyn = arb::mechanism_desc("exp2syn");
    expsyn.set("tau1", 0.5);
    expsyn.set("tau2", 1.5);
    expsyn.set("e", 0);

    cell.add_synapse({0, 0.5}, expsyn);
    return cell;
}

