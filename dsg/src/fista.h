#pragma once

#include "graph.h" //Include the graph class
#include "stats.h" //Include the stats class

class FISTA
{
private:
    double best_density = 0.0;

public:
    /**
     * @brief Perform the FISTA algorithm on the loaded graph
     */
    void
    run(const Graph &G, Stats &stats, int T = 100);

    double get_density() const
    {
        return best_density;
    }

    /**
     * @brief Extract the density of the densest subgraph from the current load vector.
     */
    double extract_density(const std::vector<int> &indices, const std::vector<std::vector<int>> &bidirectional_list, const int n, const int m, const double *b);
};
