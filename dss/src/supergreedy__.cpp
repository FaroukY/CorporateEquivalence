#include <algorithm>
#include <cmath>
#include <cstdio>
#include <execution>
#include <numeric>
#include <tbb/tbb.h>
#include <vector>

#include "graph.h"
#include "logger.h"
#include "stats.h"
#include "supergreedy.h"

/**
 * Implementation of SuperGreedy++ for the generalized p-mean densest subgraph problem.
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
 */

double get_numerator(const std::vector<double> &degrees, double p)
{
    return std::transform_reduce(
        std::execution::par, degrees.begin(), degrees.end(), 0.0,
        std::plus<double>(), [p](double x) { return std::pow(x, p); });
}

void SUPERGREEDY::run(const Graph &G, Stats &stats, const std::vector<double> y, double p, int T)
{
    const auto adj_ = G.get_bidirectional_list();
    int n = G.get_num_nodes();

    // Initialize degrees
    std::vector<double> degrees_(n, 0);
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        degrees_[i] = adj_.at(i).size();
    }

    // Initialize power table for p
    std::vector<double> power_table(n);
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        power_table[i] = std::pow(i, p);
    }

    const auto numerator_ = get_numerator(degrees_, p);
    double best_density = std::pow(static_cast<double>(numerator_) / n, 1.0 / p);

    std::vector<double> loads(n, 0);
    std::vector<bool> deleted(n, false);
    std::vector<double> D(n, 0.0);
    for (int t = 0; t < T; t++)
    {
        stats.start_timer();
        int denominator = n;
        auto numerator = numerator_;
        auto prev_numerator = numerator_;
        auto adj = adj_;
        auto degrees = degrees_;
        std::fill(deleted.begin(), deleted.end(), false);
        std::fill(D.begin(), D.end(), 0.0);

// Initialize marginals + loads
#pragma omp parallel for
        for (int i = 0; i < n; i++)
        {
            D[i] = loads[i] + std::pow(degrees[i], p);
            for (const auto &nb : adj[i])
            {
                D[i] += std::pow(degrees[nb], p) - std::pow(degrees[nb] - 1, p);
            }
        }

        // Start greedy peeling
        int next = 1;
        denominator -= 1;
        while (next <= n)
        {
            // Find min node i
            const auto i = std::min_element(
                               std::execution::par, D.begin(), D.end()) -
                           D.begin();

            // Remove i from the graph
            const auto Ni = adj[i];     // Neighbors of i
            const auto di = degrees[i]; // Degree of i

            loads[i] = D[i];
            degrees[i] = 0;
            deleted[i] = true;
            D[i] = MAXFLOAT; // So it won't be selected again

            if (di > 0)
            {
                // The numerator lost i's direct contribution.
                numerator -= power_table[di];
                // Every neighbor of i will lose the change induced by i
                const auto dichange = power_table[di] - power_table[di - 1];

// Now account for Di change in neighbors
#pragma omp parallel for
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
#pragma omp atomic
                        D[nn] += ndbchange;
                    }
                    degrees[nb]--;
#pragma omp atomic
                    // Remove 1 from degree for each neighbor of i in numerator
                    numerator -= dbchange;
                }
            }

            best_density = std::max(best_density,
                                    std::pow(static_cast<double>(numerator) / std::max(denominator, 1), 1.0 / p));
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
        log(" ");

        stats.push(t, best_density, load_vector_norm);
    }
}