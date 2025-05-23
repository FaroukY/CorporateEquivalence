#pragma once

#include "graph.h" //Include the graph class
#include "stats.h" //Include the stats class

class FRANKWOLFE
{
public:
    // Perform the FISTA algorithm on the loaded graph
    void run(const Graph &G, Stats &stats, int source, int sink, int T = 100);

    /**
     * @brief Extract the density of the densest subgraph from the current load vector.
     */
    double extract_density(int source, int sink, std::vector<int> in, std::vector<int> out, const std::vector<int> &indices, const std::vector<std::vector<std::pair<int, double>>> &bidirectional_list, const int n) const;
};
