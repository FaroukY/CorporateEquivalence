#include "frankwolfe.h"
#include "mnp.h"
#include "parser.h"
#include "stats.h"
#include "supergreedy.h"
#include <fstream>
#include <iostream>
#include <string>

std::vector<double> parse_vector_from_file(const std::string &filename, int index)
{
    std::ifstream file(filename);
    std::string line;
    int count = 0;

    // Skip the first two lines
    std::getline(file, line);
    std::getline(file, line);

    while (std::getline(file, line))
    {
        if (count == index)
        {
            std::size_t colon_pos = line.find(':');
            if (colon_pos != std::string::npos)
            {
                std::string vector_part = line.substr(colon_pos + 1);
                std::istringstream vec_stream(vector_part);
                std::vector<double> result;
                double val;
                while (vec_stream >> val)
                {
                    result.push_back(val);
                }
                return result;
            }
        }
        count++;
    }

    throw std::runtime_error("Index out of range or file format error");
}

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
    double p = 1;         // default
    int contra_query_idx = -1;
    if (algorithm == "mnp" || algorithm == "supergreedy" || algorithm == "frankwolfe")
    {
        if (argc >= 5)
        {
            iterations = std::stoi(argv[4]);
        }
        if (argc >= 6)
        {
            p = std::stod(argv[5]);
        }
        if (argc >= 7)
        {
            contra_query_idx = std::stoi(argv[6]);
        }
    }
    std::cout << "Algorithm: " << algorithm << ", iterations: " << iterations << ", p: " << p << std::endl;

    // Step 1: Read and parse the graph
    Parser parser;
    Graph graph = parser.read_graph_from_structured_file(graph_path); // assume it returns Graph
    std::cout << "Graph loaded. n = " << graph.get_num_nodes()
              << ", m = " << graph.get_num_edges() << std::endl;

    // Use query of 0 by default, this is basic p-mean DSG
    std::vector contra_vector(graph.get_num_nodes(), 0.0);
    bool is_contra = contra_query_idx > -1;
    if (is_contra)
    {
        // Fetch the query if one is provided at path <datasetname>_perturbed.txt
        std::string contra_path = graph_path;
        std::size_t dot_pos = contra_path.find_last_of('.');
        if (dot_pos != std::string::npos)
        {
            contra_path = contra_path.substr(0, dot_pos) + "_perturbed.txt";
        }
        std::cout << "Attempting to parse vector from file: " << contra_path << " with idx " << contra_query_idx << std::endl;

        try
        {
            contra_vector = parse_vector_from_file(contra_path, contra_query_idx);
        } catch (const std::exception &e)
        {
            std::cerr << "Error parsing vector from file: " << e.what() << std::endl;
            return 1;
        }
    }

    auto stats = Stats(iterations);
    if (algorithm == "supergreedy")
    {
        SUPERGREEDY sg;
        sg.run(graph, stats, contra_vector, p, iterations, is_contra);
    } else if (algorithm == "mnp")
    {
        MNP mnp;
        mnp.run(graph, stats, contra_vector, p, iterations, is_contra);
    } else if (algorithm == "frankwolfe")
    {
        FRANKWOLFE fw;
        fw.run(graph, stats, contra_vector, p, iterations, is_contra);
    } else
    {
        std::cerr << "Unknown algorithm: " << algorithm << std::endl;
        return 1;
    }

    stats.dump(output_path);
    return 0;
}
