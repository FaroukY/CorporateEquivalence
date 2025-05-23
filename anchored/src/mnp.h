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
    void run(const Graph &G, Stats &stats, const std::vector<int> &penalty, int T = 100);

    std::tuple<Eigen::VectorXd, double> min_norm_alg(const std::vector<int> &penalty, Stats &stats, int numerator, const std::vector<int> &degrees, const std::vector<std::vector<int>> &adj,
                                                     int n, int m, int iterations);

    /**
     * @brief Extract the density of the densest subgraph from the current load vector.
     */
    double extract_density(const std::vector<int> &penalty, const std::vector<int> &indices, int numerator, std::vector<int> degrees, const std::vector<std::vector<int>> &adj, const int n, const int m) const;
};