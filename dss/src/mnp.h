#pragma once

#include <Eigen/Cholesky>
#include <Eigen/Dense>

#include "graph.h" //Include the graph class
#include "stats.h" //Include the stats class

class MNP
{
private:
    double best_density_ = 0.0;
public:
    // Perform the MNP algorithm on the loaded graph
    // This is an implementation of Algorithm 1 from
    // "Provable Submodular Minimization using Wolfe’s Algorithm"
    // by Chakrabarty et al. (NIPS 2014)
    void run(const Graph &G, Stats &stats, const std::vector<double> &y, double p, int T = 100, bool is_contra = false);

    std::tuple<Eigen::VectorXd, double> min_norm_alg(Stats &stats, const std::vector<double> &y, double p, double numerator, const std::vector<double> &degrees, const std::vector<std::vector<int>> &adj,
                                                     int n, int m, int iterations, bool is_contra);

    /**
     * @brief Extract the density of the densest subgraph from the current load vector.
     */
    double extract_density(const std::vector<int> &indices, const std::vector<double> &y, double p,
                           double numerator,
                           std::vector<double> degrees, const std::vector<std::vector<int>> &adj, const int n, const int m, bool is_contra = false) const;
};