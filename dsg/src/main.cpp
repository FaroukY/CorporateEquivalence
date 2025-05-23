#include "fista.h"
#include "frankwolfe.h"
#include "greedypp.h"
#include "incremental.h"
#include "mnp.h"
#include "parser.h"
#include "pushrelabel.h"
#include "pushrelabel_contra.h"
#include "rcdm_permutation.h"
#include "stats.h"
#include <iostream>
#include <string>

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        std::cerr << "Usage: " << argv[0] << " <graph_file_path> <algorithm_name> [iterations]" << std::endl;
        return 1;
    }

    std::string graph_path = argv[1];
    std::string algorithm = argv[2];
    std::string output_path = argv[3];

    std::cout << "output path: " << output_path << std::endl;

    // Optional: handle iteration count
    int iterations = 100; // default
    if (algorithm == "mnp" || algorithm == "fista" || algorithm == "frankwolfe" || algorithm == "rcdm" || algorithm == "greedypp")
    {
        if (argc >= 5)
        {
            iterations = std::stoi(argv[4]);
        }
    }
    std::cout << "Algorithm: " << algorithm << ", iterations: " << iterations << std::endl;

    // Step 1: Read and parse the graph
    Parser parser;
    Graph graph = parser.read_graph_from_structured_file(graph_path); // assume it returns Graph
    std::cout << "Graph loaded. n = " << graph.get_num_nodes()
              << ", m = " << graph.get_num_edges() << std::endl;

    // Step 2: Dispatch to the algorithm
    auto stats = Stats(iterations);
    if (algorithm == "mnp")
    {
        MNP mnp;
        mnp.run(graph, stats, iterations);
    } else if (algorithm == "fista")
    {
        FISTA fista;
        fista.run(graph, stats, iterations);
    } else if (algorithm == "frankwolfe")
    {
        FRANKWOLFE frankwolfe;
        frankwolfe.run(graph, stats, iterations);
    } else if (algorithm == "rcdm")
    {
        runRCDM(graph, iterations, stats);
    } else if (algorithm == "incremental")
    {
        // Redirect input file to stdin
        if (freopen(graph_path.c_str(), "r", stdin) == nullptr)
        {
            std::cerr << "Error: could not open " << graph_path << " for reading." << std::endl;
            return 1;
        }
        // Fake argv so incremental gets "program name" and no arguments
        char *argv_fake[1] = {(char *)"incremental"};
        run_incremental(1, argv_fake, stats);
    } else if (algorithm == "pushrelabel")
    {
        // Redirect input file to stdin
        if (freopen(graph_path.c_str(), "r", stdin) == nullptr)
        {
            std::cerr << "Error: could not open " << graph_path << " for reading." << std::endl;
            return 1;
        }
        // pick your accuracy however you like:
        int accuracy = 100000;

        // build a little argv array with program name + accuracy
        std::string acc_str = std::to_string(accuracy);
        char *argv_fake[2];
        argv_fake[0] = (char *)"pushrelabel";
        argv_fake[1] = (char *)acc_str.c_str();

        // now argc_fake is 2, and argv_fake[1] is "100"
        run_pushrelabel(2, argv_fake, stats, graph);
    } else if (algorithm == "pushrelabel_contra")
    {
        // Redirect input file to stdin
        if (freopen(graph_path.c_str(), "r", stdin) == nullptr)
        {
            std::cerr << "Error: could not open " << graph_path << " for reading." << std::endl;
            return 1;
        }
        // pick your accuracy however you like:
        int accuracy = 100000;

        // build a little argv array with program name + accuracy
        std::string acc_str = std::to_string(accuracy);
        char *argv_fake[2];
        argv_fake[0] = (char *)"pushrelabel";
        argv_fake[1] = (char *)acc_str.c_str();

        // now argc_fake is 2, and argv_fake[1] is "100"

        std::string outfile_name;
        if (graph_path.size() >= 4 && graph_path.substr(graph_path.size() - 4) == ".txt")
        {
            outfile_name = graph_path.substr(0, graph_path.size() - 4) + "_contra.txt";
        }
        run_pushrelabel_contra(2, argv_fake, stats, outfile_name);
    } else if (algorithm == "greedypp")
    {
        GreedyPeelingPP greedypp;
        greedypp.run(graph, stats, iterations); // optional: you can limit iterations to 1, but keeping param for now
    } else
    {
        std::cerr << "Unknown algorithm: " << algorithm << std::endl;
        return 1;
    }

    stats.dump(output_path);
    return 0;
}
