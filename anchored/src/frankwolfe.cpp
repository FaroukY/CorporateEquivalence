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
    double get_numerator(const std::vector<int> &degrees, const std::vector<int> &penalty)
    {
        const auto degree_sum = std::accumulate(
            degrees.begin(), degrees.end(), 0.0);
        // Initially all elements are included.
        const auto penalty_sum = std::accumulate(
            penalty.begin(), penalty.end(), 0.0);
        return degree_sum - penalty_sum; // 2e[G] - Vol(G n R')
    }

    double norm_load_vector(const std::vector<double> &b, const int n)
    {
        return std::inner_product(b.begin(), b.end(), b.begin(), 0.0);
    }

    void LO(const std::vector<int> &penalty, std::vector<int> &indices,
            const std::vector<std::vector<int>> &adj,
            std::vector<double> &b, std::vector<double> &q,
            std::vector<bool> &deleted, std::vector<int> degrees)
    {
        int n = static_cast<int>(adj.size());
        // Reset deletion flags and direction vector
        std::fill(std::execution::par, deleted.begin(), deleted.end(), false);
        std::fill(std::execution::par, q.begin(), q.end(), 0.0);

        // Sort nodes in ascending order of current b-values to determine the extreme point
        std::sort(std::execution::par, indices.begin(), indices.end(), [&b](int aa, int bb) { return b[aa] < b[bb]; });
        // Set q(vi) = f(Si) - f(Si-1)
        auto denominator = n;
        for (int i = 0; i < indices.size(); i++)
        {
            auto idx = indices[i];
            auto numerator_delta = 2 * degrees[idx] - penalty[idx];
            for (auto nbi = 0; nbi < adj[idx].size(); nbi++)
            {
                const auto nb = adj[idx][nbi];
                if (deleted[nb])
                    continue; // Skip deleted nodes
                degrees[nb]--;
            }
            q[idx] = numerator_delta;
            denominator--;
            deleted[idx] = true;
            degrees[idx] = 0;
        }
    }
} // namespace

// Run the Frank-Wolfe algorithm for the fractional densest subgraph problem
void FRANKWOLFE::run(const Graph &graph, Stats &stats, const std::vector<int> &penalty, int T)
{
    // Number of vertices and edges in the graph
    int n = graph.get_num_nodes();
    int m = graph.get_num_edges();

    // Adjacency list: vector of (neighbor, edge_id) pairs; here we ignore the edge_id
    const auto adj_ = graph.get_bidirectional_list();

    std::vector<int> degrees_(n, 0);

#pragma omp parallel for
    for (int i = 0; i < n; ++i)
    {
        degrees_[i] = adj_.at(i).size();
    }
    const auto numerator_ = get_numerator(degrees_, penalty);

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
        LO(penalty, indices, adj_, b, q, deleted, degrees_);

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
        const auto best_density = extract_density(penalty, indices, numerator_, degrees_, adj_, n, m);

        log("Iteration ", t, ": ");
        log("Norm value: ", load_vector_norm);
        log("Best density: ", best_density);

        stats.push(t, best_density, load_vector_norm);
    }

    // After fractional solution, perform "fractional peeling" to extract an integral subgraph
    // Sort nodes by final b-value for peeling order
    std::sort(indices.begin(), indices.end(), [&](int i, int j) { return b[i] < b[j]; });
    const auto best_density = extract_density(penalty, indices, numerator_, degrees_, adj_, n, m);

    // Output the best density found
    log("Final stats: ");
    log("Norm value: ", norm_load_vector(b, n));
    log("Final density: ", best_density);
}

double FRANKWOLFE::extract_density(const std::vector<int> &penalty, const std::vector<int> &indices, int numerator, std::vector<int> degrees, const std::vector<std::vector<int>> &adj, const int n, const int m) const
{
    double highest_density = static_cast<double>(numerator) / n;
    std::vector<bool> deleted(n, false);

    auto denominator = n;
    for (int i = 0; i < indices.size(); i++)
    {
        auto idx = indices[i];
        auto numerator_delta = -2 * degrees[idx] + penalty[idx];
        for (auto nbi = 0; nbi < adj[idx].size(); nbi++)
        {
            const auto nb = adj[idx][nbi];
            if (deleted[nb])
                continue; // Skip deleted nodes
            degrees[nb]--;
        }
        numerator = numerator + numerator_delta;
        denominator--;
        deleted[idx] = true;
        degrees[idx] = 0;

        highest_density = std::max(highest_density, static_cast<double>(numerator) / denominator);
    }

    return highest_density;
}