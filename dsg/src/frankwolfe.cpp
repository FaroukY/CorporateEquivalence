#include <algorithm>
#include <cmath>
#include <cstring>
#include <execution> // Requires C++17 or later
#include <numeric>
#include <tbb/tbb.h>
#include <utility>

#include "frankwolfe.h"
#include "graph.h"
#include "logger.h"

namespace {
    double norm_load_vector(const std::vector<double> &b, const int n)
    {
        return std::inner_product(b.begin(), b.end(), b.begin(), 0.0);
    }

    void LO(std::vector<int> &indices,
            const std::vector<std::vector<std::pair<int, int>>> &adjacency_list,
            std::vector<double> &b, std::vector<double> &y,
            std::vector<bool> &deleted)
    {
        // Reset deletion flags and direction vector
        std::fill(std::execution::par, deleted.begin(), deleted.end(), false);
        std::fill(std::execution::par, y.begin(), y.end(), 0.0);

        // Sort nodes in ascending order of current b-values to determine the extreme point
        std::sort(std::execution::par, indices.begin(), indices.end(), [&b](int aa, int bb) { return b[aa] < b[bb]; });

        // For each node in sorted order, count edges to undeleted neighbors to build y, then delete the node.
        // We can't parallelize this computation because deleted nodes can't be accessed in parallel.
        for (int idx : indices)
        {
            for (const auto &[neighbor, _] : adjacency_list[idx])
            {
                if (!deleted[neighbor])
                    y[idx] += 1; // accumulate contribution from each still-active neighbor
            }
            deleted[idx] = true; // mark node as removed for subsequent counts
        }
    }
} // namespace

// Run the Frank-Wolfe algorithm for the fractional densest subgraph problem
void FRANKWOLFE::run(const Graph &graph, Stats &stats, int T)
{
    // Number of vertices and edges in the graph
    int n = graph.get_num_nodes();
    int m = graph.get_num_edges();

    // Adjacency list: vector of (neighbor, edge_id) pairs; here we ignore the edge_id
    std::vector<std::vector<std::pair<int, int>>> adjacency_list = graph.get_adjacency_list();

    // b[u]: current fractional selection score for each node (initialized to half its degree)
    std::vector<double> b(n, 0.0);
#pragma omp parallel for
    for (int u = 0; u < n; ++u)
    {
        b[u] = 0.5 * adjacency_list[u].size(); // start at mid-point between in/out for each node
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
        // Set y to the direction vector from the oracle
        LO(indices, adjacency_list, b, y, deleted);

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
        const auto best_density = extract_density(indices, graph.get_bidirectional_list(), n, m);

        log("Iteration ", t, ": ");
        log("Norm value: ", load_vector_norm);
        log("Best density: ", best_density);

        stats.push(t, best_density, load_vector_norm);
    }

    // After fractional solution, perform "fractional peeling" to extract an integral subgraph
    // Sort nodes by final b-value for peeling order
    std::sort(indices.begin(), indices.end(), [&](int i, int j) { return b[i] < b[j]; });
    const auto best_density = extract_density(indices, graph.get_bidirectional_list(), n, m);

    // Output the best density found
    log("Final stats: ");
    log("Norm value: ", norm_load_vector(b, n));
    log("Final density: ", best_density);
}

double FRANKWOLFE::extract_density(const std::vector<int> &indices, const std::vector<std::vector<int>> &bidirectional_list, const int n, const int m) const
{
    int num_nodes = n;
    int num_edges = m;
    bool *deleted = new bool[n];
    std::fill(deleted, deleted + n, false);

    // Track the best density seen during peeling
    double best_density = static_cast<double>(m) / n;

    // Remove nodes one by one in increasing order of b, updating edge count
    for (int i = 0; i < n; ++i)
    {
        int u = indices[i];
        num_nodes--; // one fewer node in the subgraph
        for (auto v : bidirectional_list[u])
        {
            if (!deleted[v])
                num_edges--; // remove edges incident to the peeled node
        }
        // Compute current density if subgraph is non-empty
        if (num_nodes > 0)
        {
            double density = static_cast<double>(num_edges) / num_nodes;
            if (density > best_density)
                best_density = density;
        }
        deleted[u] = true; // mark node as peeled
    }

    return best_density;
}