#pragma once

#include "graph.h" //Include the graph class
#include "stats.h" //Include the stats class

class FRANKWOLFE
{
private:
    double best_density_ = 0.0;
public:
    // Perform the FISTA algorithm on the loaded graph
    void run(const Graph &G, Stats &stats, const std::vector<double> &y, double p, int T = 100, bool is_contra = false);

    /**
     * @brief Extract the density of the densest subgraph from the current load vector.
     */
    double extract_density(const std::vector<int> &indices, const std::vector<double> &y, double p,
                           double numerator,
                           std::vector<double> degrees, const std::vector<std::vector<int>> &bidirectional_list, const int n, bool is_contra) const;
};
