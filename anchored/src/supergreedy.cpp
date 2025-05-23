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

double get_numerator(const std::vector<int> &degrees, const std::vector<int> &penalty)
{
    const auto degree_sum = std::accumulate(
        degrees.begin(), degrees.end(), 0.0);
    // Initially all elements are included.
    const auto penalty_sum = std::accumulate(
        penalty.begin(), penalty.end(), 0.0);
    return degree_sum - penalty_sum; // 2e[G] - Vol(G n R')
}

void SUPERGREEDY::run(const Graph &G, Stats &stats, const std::vector<int> &penalty, int T)
{
    const auto adj = G.get_bidirectional_list();
    int n = G.get_num_nodes();

    // Initialize degrees
    std::vector<int> degrees_(n, 0);
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        degrees_[i] = adj[i].size();
    }

    const auto numerator_ = get_numerator(degrees_, penalty);
    double best_density = static_cast<double>(numerator_) / n;

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
            // Removal of node i removes its penalty but also
            // removes its own degree contribution, as well as 1 from each
            // adjacent neighbour.
            D[i] = loads[i] - penalty[i] + 2 * degrees[i];
        }

        // Initialize heap
        heap.clear();
        for (int i = 0; i < n; i++)
        {
            handles[i] = heap.push({i, D[i]});
        }

        // Start greedy peeling
        while (!heap.empty())
        {
            auto min = heap.top();
            const auto i = min.key;
            const auto Di = min.priority;
            int diff = Di - loads[i];
            loads[i] = Di;              // Store the argmin marginal + load
            deleted[i] = true;          // Mark i as deleted
            heap.pop();

            // Remove i's penalty from the numerator
            numerator -= diff;
            denominator -= 1;
            if (denominator)
                best_density = std::max(best_density, static_cast<double>(numerator) / denominator);

            for (auto nbi = 0; nbi < adj[i].size(); nbi++)
            {
                const auto nb = adj[i][nbi];
                if (deleted[nb])
                    continue; // Skip deleted nodes

                // nb has lost a degree, so we remove 1.
                // AKA: removal of nb now incurs a cost that's less by 2 * (1 edge)
                D[nb] = D[nb] - 2;
                heap.decrease(handles[nb], {nb, D[nb]});
            }
        }
        stats.pause_timer();
        const auto load_vector_norm = std::transform_reduce(
            loads.begin(), loads.end(), 0.0,
            std::plus<double>(), [t](double x) { return std::pow(x / (t + 1), 2); });

        std::cout << "Best_density: " << best_density << std::endl;
        log("numerator: ", numerator);
        log("Iteration ", t, ": ");
        log("Best density: ", best_density);
        log("Load norm: ", load_vector_norm);
        log("");

        stats.push(t, best_density, load_vector_norm);
    }
}