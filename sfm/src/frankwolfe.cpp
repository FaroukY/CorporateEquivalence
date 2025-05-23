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
    double get_numerator(const std::vector<double> &degrees, double p)
    {
        return std::transform_reduce(
            std::execution::par, degrees.begin(), degrees.end(), 0.0,
            std::plus<double>(), [p](double x) { return std::pow(x, p); });
    }

    double norm_load_vector(const std::vector<double> &b, const int n)
    {
        return std::inner_product(b.begin(), b.end(), b.begin(), 0.0);
    }

    void LO(int source, int sink, std::vector<int> &indices,
            const std::vector<std::vector<std::pair<int, double>>> &adj,
            std::vector<double> &b, std::vector<double> &y,
            std::vector<bool> &deleted, std::vector<int> in, std::vector<int> out)
    {
        int n = static_cast<int>(adj.size());
        // Reset deletion flags and direction vector
        std::fill(std::execution::par, deleted.begin(), deleted.end(), false);
        std::fill(std::execution::par, y.begin(), y.end(), 0.0);

        // Sort nodes in ascending order of current b-values to determine the extreme point
        std::sort(std::execution::par, indices.begin(), indices.end(), [&b](int aa, int bb) { return b[aa] < b[bb]; });

        deleted[source] = true;
        deleted[sink] = true;

        // Set q(vi) = f(Si) - f(Si-1)
        for (int i = 0; i < indices.size(); i++)
        {
            if (deleted[indices[i]])
            {
                continue;
            }

            auto idx = indices[i];
            // Peel
            const auto delta = in[idx] - out[idx];
            y[idx] = delta;
            deleted[idx] = true;
            // Update neighbors
            for (const auto &[nbi, w] : adj[idx])
            {
                if (deleted[nbi])
                {
                    continue;
                }
                // Update the score. The outgoing will increase by w, the inner will decrease by w.
                // So overall the score increases by 2 * w
                out[nbi] += w;
                in[nbi] -= w;
            }
        }
    }
} // namespace

// Run the Frank-Wolfe algorithm for the fractional densest subgraph problem
void FRANKWOLFE::run(const Graph &graph, Stats &stats, int source, int sink, int T)
{
    // Number of vertices and edges in the graph
    int n = graph.get_num_nodes();
    int m = graph.get_num_edges();

    // Adjacency list: vector of (neighbor, edge_id) pairs; here we ignore the edge_id
    const auto adj_ = graph.get_bidirectional_list();

    std::vector<int> out_(n, 0);
    std::vector<int> in_(n, 0);

    // Initialize out and in degrees for all base nodes. Initially
    // out is t and in is everything else.
    for (int i = 0; i < n; ++i)
    {
        for (const auto &[nbi, w] : adj_[i])
        {
            if (nbi == sink)
            {
                out_[i] += w;
            } else
            {
                in_[i] += w;
            }
        }
    }

    // b[u]: current fractional selection score for each node
    std::vector<double> b(n, 0.0);
#pragma omp parallel for
    for (int u = 0; u < n; ++u)
    {
        b[u] = 0.5 * adj_[u].size(); // start at mid-point between in/out for each node
    }

    // y[u]: temporary direction vector from the linear minimization oracle
    std::vector<double> y(n, 0.0);
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
        LO(source, sink, indices, adj_, b, y, deleted, in_, out_);

// --- Update step ---
// Move b toward the oracle direction y by step size gamma
#pragma omp parallel for
        for (int i = 0; i < n; ++i)
        {
            b[i] = (1 - gamma) * b[i] + gamma * y[i];
        }

        // Calculate and store stats
        stats.pause_timer();
        std::sort(std::execution::par, indices.begin(), indices.end(), [&b](int i, int j) { return b[i] < b[j]; });

        const auto load_vector_norm = norm_load_vector(b, n);
        const auto best_density = extract_density(source, sink, in_, out_, indices, graph.get_bidirectional_list(), n);

        log("Iteration ", t, ": ");
        log("Norm value: ", load_vector_norm);
        log("Best density: ", best_density);

        stats.push(t, best_density, load_vector_norm);
    }

    // After fractional solution, perform "fractional peeling" to extract an integral subgraph
    // Sort nodes by final b-value for peeling order
    std::sort(indices.begin(), indices.end(), [&](int i, int j) { return b[i] < b[j]; });
    const auto best_density = extract_density(source, sink, in_, out_, indices, graph.get_bidirectional_list(), n);

    // Output the best density found
    log("Final stats: ");
    log("Norm value: ", norm_load_vector(b, n));
    log("Final density: ", best_density);
}

double FRANKWOLFE::extract_density(int source, int sink, std::vector<int> in, std::vector<int> out, const std::vector<int> &indices, const std::vector<std::vector<std::pair<int, double>>> &adj, const int n) const
{
    double best_cut = -in[sink];
    double current_cut = best_cut;
    std::vector<bool> deleted(n, false);
    deleted[sink] = true;
    deleted[source] = true;

    for (int i = 0; i < indices.size(); i++)
    {
        auto idx = indices[i];
        if (deleted[idx])
        {
            continue;
        }
        const auto delta = in[idx] - out[idx];
        current_cut -= delta;
        deleted[idx] = true;
        // Update neighbors
        for (const auto &[nbi, w] : adj[idx])
        {
            if (deleted[nbi])
            {
                continue;
            }
            // Update the score. The outgoing will increase by w, the inner will decrease by w.
            // So overall the score increases by 2 * w
            out[nbi] += w;
            in[nbi] -= w;
        }
        best_cut = std::max(best_cut, current_cut);
    }
    return best_cut;
}