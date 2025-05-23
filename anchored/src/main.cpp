#include "frankwolfe.h"
#include "mnp.h"
#include "parser.h"
#include "stats.h"
#include "supergreedy.h"
#include "flow.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_set>

std::vector<int> parse_vector(const std::string &filename)
{
    std::vector<int> numbers;
    std::ifstream infile(filename);

    if (!infile.is_open())
    {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return numbers;
    }

    std::string line;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        int number;
        if (iss >> number)
        {
            numbers.push_back(number);
        }
    }

    return numbers;
}

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        std::cerr << "Usage: " << argv[0] << " <graph_file_path> <algorithm_name> [iterations]" << std::endl;
        return 1;
    }

    std::string penalty_path = argv[1]; // We get the graph path from this
    std::string algorithm = argv[2];
    std::string output_path = argv[3];

    std::cout << "output path: " << output_path << std::endl;
    
    // Optional: handle iteration count
    int iterations = 100; // default
    int penalty_idx = -1;
    if (algorithm == "mnp" || algorithm == "supergreedy" || algorithm == "frankwolfe")
    {
        if (argc >= 5)
        {
            iterations = std::stoi(argv[4]);
        }
    }
    std::cout << "Algorithm: " << algorithm << ", iterations: " << iterations << std::endl;

    // Step 1: Read and parse the graph
    Parser parser;

    // Graph path is penalty up to _R
    std::string graph_path = penalty_path.substr(0, penalty_path.find_first_of('_')) + ".txt";
    std::cout << "Graph path: " << graph_path << std::endl;

    Graph graph = parser.read_graph_from_structured_file(graph_path); // assume it returns Graph
    std::cout << "Graph loaded. n = " << graph.get_num_nodes()
              << ", m = " << graph.get_num_edges() << std::endl;

    std::vector<int> R_vec;
    if (penalty_path.find(".txt") != std::string::npos)
    {
        R_vec = parse_vector(penalty_path);
        std::cout << "Penalty loaded. Size: " << R_vec.size() << std::endl;
    } else
    {
        std::cerr << "Invalid penalty file format. Expected .txt" << std::endl;
        return 1;
    }
    std::unordered_set<int> R_set(R_vec.begin(), R_vec.end());

    // The cost of each node NOT in the R_set is its degree in G. Nodes
    // in the R_set have a cost of 0.
    std::cout << "R_set size: " << R_set.size() << std::endl;
    std::vector<int> penalty(graph.get_num_nodes(), 0);
    auto adj = graph.get_adjacency_list();
    for (int i = 0; i < graph.get_num_nodes(); ++i)
    {
        if (R_set.find(i) == R_set.end())
        {
            penalty[i] = adj[i].size();
        }
    }

    auto stats = Stats(iterations);
    if (algorithm == "supergreedy")
    {
        SUPERGREEDY sg;
        sg.run(graph, stats, penalty, iterations);
    } else if (algorithm == "mnp")
    {
        MNP mnp;
        mnp.run(graph, stats, penalty, iterations);
    } else if (algorithm == "frankwolfe")
    {
        FRANKWOLFE fw;
        fw.run(graph, stats, penalty, iterations);
    } else if (algorithm == "flow"){
        run_flow(graph, stats, penalty);
    } 
    else
    {
        std::cerr << "Unknown algorithm: " << algorithm << std::endl;
        return 1;
    }

    stats.dump(output_path);
    return 0;
}
