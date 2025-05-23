#include <algorithm>
#include <cmath>
#include <cstring>
#include <execution> // Requires C++17 or later
#include <numeric>
#include <tbb/tbb.h>

#include "frankwolfe.h"
#include "graph.h"
#include "logger.h"

namespace {
    double get_numerator(const std::vector<double> &y, const std::vector<double> &degrees, double p)
    {
        const auto degree_sum = std::transform_reduce(
            std::execution::par, degrees.begin(), degrees.end(), 0.0,
            std::plus<double>(), [p](double x) { return std::pow(x, p); });
        const auto y_sum = std::accumulate(
            y.begin(), y.end(), 0.0);
        return degree_sum - 2 * y_sum;
    }

    double norm_load_vector(const std::vector<double> &b, const int n)
    {
        return std::inner_product(b.begin(), b.end(), b.begin(), 0.0);
    }

    void LO(const std::vector<double> &y, std::vector<int> &indices,
            const std::vector<std::vector<int>> &adj,
            std::vector<double> &b, std::vector<double> &q,
            std::vector<bool> &deleted, double p, double numerator, std::vector<double> degrees)
    {
        int n = static_cast<int>(adj.size());
        // Reset deletion flags and direction vector
        std::fill(std::execution::par, deleted.begin(), deleted.end(), false);
        std::fill(std::execution::par, q.begin(), q.end(), 0.0);

        // Sort nodes in ascending order of current b-values to determine the extreme point
        std::sort(std::execution::par, indices.begin(), indices.end(), [&b](int aa, int bb) { return b[aa] < b[bb]; });
        // Set q(vi) = f(Si) - f(Si-1)
        for (int i = 0; i < indices.size(); i++)
        {
            auto idx = indices[i];
            // log("Peeling node ", idx, " with degree ", degrees[idx], " and numerator ", numerator);
            auto numerator_delta = std::pow(degrees[idx], p) - 2 * y[idx];
            for (auto nbi = 0; nbi < adj[idx].size(); nbi++)
            {
                const auto nb = adj[idx][nbi];
                if (deleted[nb])
                    continue; // Skip deleted nodes

                if (degrees[nb] > 0)
                    numerator_delta += std::pow(degrees[nb], p) - std::pow(degrees[nb] - 1, p);
                degrees[nb]--;
            }
            // log("numerator ", new_numerator, " denominator ", denominator, " numerator_delta ", numerator_delta);
            q[idx] = numerator_delta;

            deleted[idx] = true;
            degrees[idx] = 0;
        }
    }
} // namespace

// Run the Frank-Wolfe algorithm for the fractional densest subgraph problem
void FRANKWOLFE::run(const Graph &graph, Stats &stats, const std::vector<double> &y, double p, int T, bool is_contra)
{
    // Number of vertices and edges in the graph
    int n = graph.get_num_nodes();
    int m = graph.get_num_edges();

    // Adjacency list: vector of (neighbor, edge_id) pairs; here we ignore the edge_id
    const auto adj_ = graph.get_bidirectional_list();

    std::vector<double> degrees_(n, 0);

#pragma omp parallel for
    for (int i = 0; i < n; ++i)
    {
        degrees_[i] = adj_.at(i).size();
    }
    const auto numerator_ = get_numerator(y, degrees_, p);

    // b[u]: current fractional selection score for each node (initialized to half its degree)
    std::vector<double> b(n, 0.0);
#pragma omp parallel for
    for (int u = 0; u < n; ++u)
    {
        b[u] = 0.5 * adj_[u].size(); // start at mid-point between in/out for each node
    }

    // y[u]: temporary direction vector from the linear minimization oracle
    std::vector<double> q(n, 0.0);
    // indices: list of all node IDs, used for sorting in the oracle
    std::vector<int> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    // deleted[u]: flag to mark nodes as removed during the oracle pass
    std::vector<bool> deleted(n, false);

    // Main Frank-Wolfe loop: T iterations
    for (int t = 0; t < T; ++t)
    {
        stats.start_timer();
        // Step size gamma = 1/(t+1)
        double gamma = 1.0 / (t + 1);

        // --- Linear Minimization Oracle ---
        // Compute new extreme direction y = LO(adj, b)
        LO(y, indices, adj_, b, q, deleted, p, numerator_, degrees_);

// --- Update step ---
// Move b toward the oracle direction y by step size gamma
#pragma omp parallel for
        for (int i = 0; i < n; ++i)
        {
            b[i] = (1 - gamma) * b[i] + gamma * q[i];
        }

        // Calculate and store stats
        stats.pause_timer();
        std::sort(std::execution::par, indices.begin(), indices.end(), [&b](int i, int j) { return b[i] < b[j]; });

        const auto load_vector_norm = norm_load_vector(b, n);
        best_density_ = std::max(best_density_, extract_density(indices, y, p, numerator_, degrees_, graph.get_bidirectional_list(), n, is_contra));

        log("Iteration ", t, ": ");
        log("Norm value: ", load_vector_norm);
        log("Best density: ", best_density_);

        stats.push(t, best_density_, load_vector_norm);
    }

    // After fractional solution, perform "fractional peeling" to extract an integral subgraph
    // Sort nodes by final b-value for peeling order
    std::sort(indices.begin(), indices.end(), [&](int i, int j) { return b[i] < b[j]; });
    best_density_ = std::max(best_density_, extract_density(indices, y, p, numerator_, degrees_, graph.get_bidirectional_list(), n, is_contra));

    // Output the best density found
    log("Final stats: ");
    log("Norm value: ", norm_load_vector(b, n));
    log("Final density: ", best_density_);
}

double FRANKWOLFE::extract_density(const std::vector<int> &indices, const std::vector<double> &y, double p,
                                   double numerator,
                                   std::vector<double> degrees, const std::vector<std::vector<int>> &adj, const int n, bool is_contra) const
{
    double highest_density = std::pow(numerator / n, 1.0 / p);
    std::vector<bool> deleted(n, false);

    auto denominator = n;
    auto denom = is_contra ? 1 : n;
    for (int i = 0; i < indices.size() - 1; i++)
    {
        auto idx = indices[i];
        auto numerator_delta = -std::pow(degrees[idx], p) + 2 * y[idx];
        for (auto nbi = 0; nbi < adj[idx].size(); nbi++)
        {
            const auto nb = adj[idx][nbi];
            if (deleted[nb])
                continue; // Skip deleted nodes

            if (degrees[nb] > 0)
                numerator_delta -= std::pow(degrees[nb], p) - std::pow(degrees[nb] - 1, p);
            degrees[nb]--;
        }
        numerator = numerator + numerator_delta;
        denominator--;
        deleted[idx] = true;
        degrees[idx] = 0;

        denom = is_contra ? 1 : std::max(denominator, 1);
        highest_density = std::max(highest_density,
                                   std::pow(numerator / denom, 1.0 / p));
    }

    return highest_density;
}