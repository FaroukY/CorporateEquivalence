#pragma once

#include <Eigen/Cholesky>
#include <Eigen/Dense>

#include "graph.h" //Include the graph class
#include "stats.h" //Include the stats class

class MNP
{
public:
    // Perform the MNP algorithm on the loaded graph
    // This is an implementation of Algorithm 1 from
    // "Provable Submodular Minimization using Wolfe’s Algorithm"
    // by Chakrabarty et al. (NIPS 2014)
    void run(const Graph &G, Stats &stats, int source, int sink, int T = 100);

    std::tuple<Eigen::VectorXd, double> min_norm_alg(int source, int sink, Stats &stats, const std::vector<std::vector<std::pair<int, double>>> &adj,
                                                     int n, int m, int iterations);

    /**
     * @brief Extract the density of the densest subgraph from the current load vector.
     */
    double extract_density(int source, int sink, const std::vector<int> &indices, std::vector<int> in,
                           std::vector<int> out, const std::vector<std::vector<std::pair<int, double>>> &adj, const int n, const int m) const;
};