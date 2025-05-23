#include "frankwolfe.h"
#include "mnp.h"
#include "parser.h"
#include "stats.h"
#include "supergreedy.h"
#include <iostream>
#include <limits>
#include <queue>
#include <regex>
#include <string>

double edmonds_karp_min_cut(const Graph &graph, int source, int sink)
{
    int n = graph.get_num_nodes();
    std::vector<std::unordered_map<int, double>> capacity(n);

    // Build capacity matrix from adjacency list
    for (int u = 0; u < n; ++u)
    {
        for (const auto &[v, w] : graph.get_bidirectional_list()[u])
        {
            capacity[u][v] += w;
        }
    }

    auto bfs = [&](std::vector<int> &parent) -> double {
        std::fill(parent.begin(), parent.end(), -1);
        parent[source] = -2;
        std::queue<std::pair<int, double>> q;
        q.push({source, std::numeric_limits<double>::max()});

        while (!q.empty())
        {
            int u = q.front().first;
            double flow = q.front().second;
            q.pop();

            for (const auto &[v, cap] : capacity[u])
            {
                if (parent[v] == -1 && cap > 0)
                {
                    parent[v] = u;
                    double new_flow = std::min(flow, cap);
                    if (v == sink)
                    {
                        return new_flow;
                    }
                    q.push({v, new_flow});
                }
            }
        }
        return 0;
    };

    double max_flow = 0;
    std::vector<int> parent(n);
    double new_flow;

    while ((new_flow = bfs(parent)) > 0)
    {
        max_flow += new_flow;
        int cur = sink;
        while (cur != source)
        {
            int prev = parent[cur];
            capacity[prev][cur] -= new_flow;
            capacity[cur][prev] += new_flow;
            cur = prev;
        }
    }

    return max_flow;
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

    // Optional: handle iteration count
    int iterations = 100; // default
    double p = 1;         // default
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
    }
    std::cout << "Algorithm: " << algorithm << ", iterations: " << iterations << ", p: " << p << std::endl;
    // Step 1: Read and parse the graph
    Parser parser;
    Graph graph = parser.read_graph_from_structured_file(graph_path); // assume it returns Graph
    std::cout << "Graph loaded. n = " << graph.get_num_nodes()
              << ", m = " << graph.get_num_edges() << std::endl;

    // DIMACS sets the source as the least node, the sink as the greatest node.
    auto s = 0;
    auto t = graph.get_num_nodes() - 1;

    // If the dataset name contains the word '.max' then s is prefixed by s{value}
    // and t is prefixed by t{value} so we need to parse it.
    std::regex pattern(R"(.*_s(\d+)_t(\d+))");
    std::smatch matches;
    if (graph_path.find(".max") != std::string::npos)
    {
        if (std::regex_search(graph_path, matches, pattern))
        {
            s = std::stoi(matches[1]);
            t = std::stoi(matches[2]);
            std::cout << "s = " << s << ", t = " << t << std::endl;
        } else
        {
            std::cerr << "Pattern not matched!" << std::endl;
        }
    }
    std::cout << "Source: " << s << ", Sink: " << t << std::endl;

    // Step 2: Dispatch to the algorithm
    auto stats = Stats(iterations);
    if (algorithm == "supergreedy")
    {
        SUPERGREEDY sg;
        sg.run(graph, stats, s, t, p, iterations);
    } else if (algorithm == "mnp")
    {
        MNP mnp;
        mnp.run(graph, stats, s, t, iterations);
    } else if (algorithm == "frankwolfe")
    {
        FRANKWOLFE fw;
        fw.run(graph, stats, s, t, iterations);
    } else if (algorithm == "exact")
    {
        stats.start_timer();
        double min_cut = edmonds_karp_min_cut(graph, s, t);
        stats.pause_timer();

        std::cout << "Exact Min Cut value: " << min_cut << std::endl;
        stats.push(0, min_cut, 0); // Just 1 stat entry
    }
    else
    {
        std::cerr << "Unknown algorithm: " << algorithm << std::endl;
        return 1;
    }

    stats.dump(output_path);
    return 0;
}
