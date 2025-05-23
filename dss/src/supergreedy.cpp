#include "stats.h"
#include <algorithm>
#include <boost/heap/fibonacci_heap.hpp>
#include <cmath>
#include <cstdio>
#include <execution>
#include <numeric>
#include <tbb/tbb.h>
#include <vector>

#include "graph.h"
#include "logger.h"
#include "supergreedy.h"

/**
 * Implementation of SuperGreedy++ for the generalized p-mean densest subgraph problem + query vector.
 * The objective is to find S that minimizes f(S) = ((1 / |S|) * sum_{i in S} d_i^p) ^ (1/p).
 * We notice that this is equivalent to maximizing f(S) = sum_{i in S} d_i^p.
 *
 * SuperGreedy++ selects argmin_v load_v + f(v | S_j - v) = f(S_j) - f(S_j - v).
 * To compute the marginal efficiently we notice that removal of a node v from S will result
 * in the loss of its degree contribution to the sum (d_v^p) and its neighbors will all lose a degree (sum_{u in N(v)} (d_u^p - (d_u - 1)^p)).
 * The marginal is then given by: -(d_v^p + sum_{u in N(v)} (d_u^p - (d_u - 1)^p)).
 *
 * We will maintain these marginals in the vector D, for each v. Notice that a neighbor of v, say w, will have its degree
 * considered in D for all its neighbors. So if we remove v, and w is a neighbor of v, then we need to update D for all
 * neighbors of w: namely take away d_w^p - (d_w - 1)^p, and replace by (d_w - 1)^p - (d_u - 2)^p
 * => -(d_w^p - (d_w - 1)^p) + (d_w - 1)^p + (d_u - 2)^p = -d_w^p + 2 * (d_w - 1)^p - (d_u - 2)^p
 *
 *
 * We augment the problem to query points in the contrapolymatroid of f(S) = E(S) = 0.5 * sum_{i in S} d_i.
 * Such points will satisfy that y(S) - f(S) >= 0 => 2 * y(S) - sum_{i in S} d_i >= 0.
 * Or that sum_{i in S} d_i - 2 * y(S) <= 0, where the first part is supermodular and the second part is modular. Hence
 * this function is supermodular and we can maximize it (or minimize the negative, submodular, function).
 * If the maximum is above 0, then the query point is a NO instance. It's a YES instance otherwise.
 */

struct Node
{
    int key;         // vertex
    double priority; // load
    bool operator<(const Node &other) const
    {
        return priority > other.priority; // min-heap, default is max-heap
    }
};

using FibHeap = boost::heap::fibonacci_heap<Node>;
using Handle = FibHeap::handle_type;

double get_numerator(const std::vector<double> &degrees, double p, const std::vector<double> y)
{
    const auto degree_sum = std::transform_reduce(
        std::execution::par, degrees.begin(), degrees.end(), 0.0,
        std::plus<double>(), [p](double x) { return std::pow(x, p); });
    // Initially all elements are included.
    const auto y_sum = std::accumulate(
        y.begin(), y.end(), 0.0);
    std::cout << "y_sum: " << y_sum << std::endl;
    return degree_sum - 2 * y_sum;
}

void SUPERGREEDY::run(const Graph &G, Stats &stats, const std::vector<double> y, double p, int T, bool is_contra)
{
    const auto adj = G.get_bidirectional_list();
    int n = G.get_num_nodes();

    // Initialize degrees
    std::vector<double> degrees_(n, 0);
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        degrees_[i] = adj[i].size();
    }

    // Initialize power table for p
    std::vector<double> power_table(n);
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        power_table[i] = std::pow(i, p);
    }

    const auto numerator_ = get_numerator(degrees_, p, y);
    auto denom = is_contra ? 1 : n;
    double best_density = std::pow(static_cast<double>(numerator_) / denom, 1.0 / p);

    std::vector<double> loads(n, 0);
    std::vector<bool> deleted(n, false);
    std::vector<double> D(n, 0.0);
    FibHeap heap;
    std::vector<Handle> handles(n);

    for (int t = 0; t < T; t++)
    {
        stats.start_timer();
        int denominator = n;
        auto numerator = numerator_;
        auto degrees = degrees_;
        std::fill(deleted.begin(), deleted.end(), false);
        std::fill(D.begin(), D.end(), 0.0);

// Initialize marginals + loads
#pragma omp parallel for
        for (int i = 0; i < n; i++)
        {
            D[i] = loads[i] - 2 * y[i] + std::pow(degrees[i], p);
            for (const auto &nb : adj[i])
            {
                D[i] += std::pow(degrees[nb], p) - std::pow(degrees[nb] - 1, p);
            }
        }

        // Initialize heap
        heap.clear();
        for (int i = 0; i < n; i++)
        {
            handles[i] = heap.push({i, D[i]});
        }

        // Start greedy peeling
        int next = 1;
        denominator--;
        while (next <= n)
        {
            auto min = heap.top();
            heap.pop();
            const auto i = min.key;
            const auto Di = min.priority;
            const auto Ni = adj[i];     // Neighbors of i
            const auto di = degrees[i]; // Degree of i
            loads[i] = Di;              // Store the argmin marginal + load
            deleted[i] = true;          // Mark i as deleted
            degrees[i] = 0;             // Set degree to 0

            // Remove the linear weight contribution for this node.
            numerator += 2 * y[i];
            if (di > 0)
            {
                // The numerator lost i's direct contribution.
                numerator -= power_table[di];
                // Every neighbor of i will lose the change induced by i
                const auto dichange = power_table[di] - power_table[di - 1];

                // Now account for Di change in neighbors
                // #pragma omp parallel for
                for (auto nbi = 0; nbi < adj[i].size(); nbi++)
                {
                    const auto nb = adj[i][nbi];
                    if (deleted[nb])
                        continue; // Skip deleted nodes

                    // Every neighbor has lost 1 degree for their own direct
                    // contribution and we remove the change induced by i.
                    const auto db = degrees[nb];
                    auto dbchange = power_table[db] - power_table[db - 1];
                    D[nb] = D[nb] - dichange - dbchange;

                    // And then account for neighbours of neighbours of di. The change
                    // for the neighbor, w, of i has changed due to reduction in w's degree.
                    const auto ndbchange = 2 * power_table[db - 1] - power_table[db - 2] - power_table[db];
                    for (auto nni = 0; nni < adj[nb].size(); nni++)
                    {
                        const auto nn = adj[nb][nni];
                        if (deleted[nn])
                            continue; // Skip deleted nodes
                                      // #pragma omp atomic
                        D[nn] += ndbchange;
                    }
                    degrees[nb]--;
                    // #pragma omp atomic
                    // Remove 1 from degree for each neighbor of i in numerator
                    numerator -= dbchange;
                    heap.decrease(handles[nb], {nb, D[nb]});
                }
            }

            if (denominator)
                denom = is_contra ? 1 : n;
                best_density = std::max(best_density,
                                        std::pow(static_cast<double>(numerator) / denom, 1.0 / p));
            denominator -= 1;
            next++;
        }
        stats.pause_timer();
        const auto load_vector_norm = std::transform_reduce(
            loads.begin(), loads.end(), 0.0,
            std::plus<double>(), [t](double x) { return std::pow(x / (t + 1), 2); });
        log("numerator: ", numerator);
        log("Iteration ", t, ": ");
        log("Best density: ", best_density);
        log("Load norm: ", load_vector_norm);
        log("");

        stats.push(t, best_density, load_vector_norm);
    }
}