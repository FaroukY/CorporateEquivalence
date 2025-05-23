#include <algorithm>
#include <cmath>
#include <cstring>
#include <execution>
#include <numeric>
#include <tbb/tbb.h>
#include <utility>

#include "fista.h"
#include "graph.h"
#include "logger.h"

namespace {
    double norm_load_vector(const double *b, const int n)
    {
        return std::inner_product(b, b + n, b, 0.0);
    }
} // namespace

// Run the FISTA algorithm for fractional densest subgraph
void FISTA::run(const Graph &graph, Stats &stats, int T)
{
    // Number of vertices and edges
    int n = graph.get_num_nodes();
    int m = graph.get_num_edges();
    log("Done reading graph");

    // Adjacency list: (neighbor, edge_index) pairs
    std::vector<std::vector<std::pair<int, int>>> adjacency_list = graph.get_adjacency_list();

    // Prepare arrays for edge-to-source mapping and reverse-edge lookup
    int *reverse_edge_idx = new int[2 * m];
    int *edge_src_indices = new int[2 * m];
#pragma omp parallel for
    for (int i = 0; i < n; ++i)
    {
        // For each outgoing edge of node i
        for (int j = 0; j < adjacency_list[i].size(); ++j)
        {
            int sister_idx = adjacency_list[i][j].second;
            // flip parity to get matching reverse edge index
            int idx = (sister_idx % 2 == 0 ? sister_idx + 1 : sister_idx - 1);
            edge_src_indices[idx] = i;          // record source node for this directed edge
            reverse_edge_idx[idx] = sister_idx; // map to its reverse
        }
    }

    // Compute step size: 0.5 divided by maximum node degree
    double max_degree = (*std::max_element(
                             std::execution::par,
                             adjacency_list.begin(), adjacency_list.end(),
                             [](auto &a, auto &b) { return a.size() < b.size(); }))
                            .size();
    double learning_rate = 0.5 / max_degree;
    log("Learning rate: ", learning_rate);
    log("Initial density: ", 1.0 * m / n);

    // Allocate FISTA variables: x = current iterate, y = momentum mix, z = gradient step
    double *x = new double[2 * m];
    double *y = new double[2 * m];
    double *z = new double[2 * m];
    double *b = new double[n]; // holds per-node sums

    // Sorted indices for tracking best density
    std::vector<int> indices(n);
    std::iota(indices.begin(), indices.end(), 0);

// Initialize x and y to 0.5 for all directed edges
#pragma omp parallel for
    for (int i = 0; i < 2 * m; ++i)
    {
        x[i] = 0.5;
        y[i] = 0.5;
    }
    // Zero out z and b arrays
    std::memset(z, 0, 2 * m * sizeof(double));

    double tk = 1.0; // momentum parameter

    // Buffers for new updates, allocated once
    double *new_xuv = new double[2 * m];
    double *new_y = new double[2 * m];

    // Main FISTA loop: T iterations
    for (int t = 0; t < T; ++t)
    {
        stats.start_timer();
        // 7a. Gradient computation: b[u] = sum of y over edges outgoing from u
        std::memset(b, 0, n * sizeof(double));
        for (int i = 0; i < 2 * m; ++i)
            b[edge_src_indices[i]] += y[i];

// 7b. Gradient step: z_{uv} = y_{uv} - 2 * lr * b_u
#pragma omp parallel for
        for (int i = 0; i < 2 * m; ++i)
            z[i] = y[i] - 2.0 * learning_rate * b[edge_src_indices[i]];

        // 7c. Update momentum scalar via FISTA formula
        double tknew = (1.0 + std::sqrt(1.0 + 4.0 * tk * tk)) / 2.0;

// 7d. Proximal / projection step to satisfy x_uv + x_vu = 1
// new_xuv[i] = clamp((z[i] - z[reverse_edge] + 1)/2, [0,1])
#pragma omp parallel for
        for (int i = 0; i < 2 * m; ++i)
            new_xuv[i] = std::clamp((z[i] - z[reverse_edge_idx[i]] + 1.0) / 2.0, 0.0, 1.0);

// FISTA momentum update: new_y = new_x + ((tk-1)/tknew)*(new_x - x_old) + (tk/tknew)*(new_x - y_old)
#pragma omp parallel for
        for (int i = 0; i < 2 * m; ++i)
            new_y[i] = new_xuv[i] + ((tk - 1.0) / tknew) * (new_xuv[i] - x[i]) + (tk / tknew) * (new_xuv[i] - y[i]);

// Copy updates back into x and y, update tk for next iteration
#pragma omp parallel for
        for (int i = 0; i < 2 * m; ++i)
        {
            x[i] = new_xuv[i];
            y[i] = new_y[i];
        }
        tk = tknew;

        // Calculate and store stats
        stats.pause_timer();
        std::memset(b, 0, n * sizeof(double));
        for (int i = 0; i < 2 * m; ++i)
            b[edge_src_indices[i]] += x[i];

        std::sort(std::execution::par, indices.begin(), indices.end(), [&](int i, int j) { return b[i] < b[j]; });

        const auto load_vector_norm = norm_load_vector(b, n);
        const auto best_density = extract_density(indices, graph.get_bidirectional_list(), n, m, b);

        // Optional monitoring: track max/min b
        double max_b = *std::max_element(std::execution::par, b, b + n);
        double min_b = *std::min_element(std::execution::par, b, b + n);

        log("Iteration ", t, ": ");
        log("Norm value: ", load_vector_norm);
        log("Best density: ", best_density);
        log("Max b: ", max_b, ", Min b: ", min_b);

        stats.push(t, best_density, load_vector_norm);
    }

    // Log last iteration
    log("Final stats: ");
    const auto load_vector_norm = norm_load_vector(b, n);
    const auto best_density = extract_density(indices, graph.get_bidirectional_list(), n, m, b);
    log("Norm value: ", load_vector_norm);
    log("Final density: ", best_density);

    // Clean up allocated memory
    delete[] reverse_edge_idx;
    delete[] edge_src_indices;
    delete[] x;
    delete[] y;
    delete[] z;
    delete[] b;
    delete[] new_xuv;
    delete[] new_y;
}

double FISTA::extract_density(const std::vector<int> &indices, const std::vector<std::vector<int>> &bidirectional_list, const int n, const int m, const double *b)
{
    // Fractional peeling: extract an integral subgraph from the fractional solution
    int num_nodes = n;
    int num_edges = m;
    bool *deleted = new bool[n];
    std::fill(deleted, deleted + n, false);

    best_density = std::max(1.0 * m / n, best_density);
    for (int i = 0; i < n; ++i)
    {
        int u = indices[i];
        num_nodes--;
        for (auto v : bidirectional_list[u])
        {
            if (!deleted[v])
                num_edges--;
        }
        if (num_nodes > 0)
        {
            double d = 1.0 * num_edges / num_nodes;
            if (d > best_density)
                best_density = d;
        }
        deleted[u] = true;
    }
    delete[] deleted;

    return best_density;
}
