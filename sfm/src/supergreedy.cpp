#include "supergreedy.h"
#include "logger.h"

#include <boost/heap/fibonacci_heap.hpp>
#include <cmath>
#include <numeric>

struct Node
{
    int key;            // vertex
    long long priority; // load (changed to long long)
    bool operator<(const Node &other) const
    {
        return priority > other.priority; // min-heap, default is max-heap
    }
};

using FibHeap = boost::heap::fibonacci_heap<Node>;
using Handle = FibHeap::handle_type;

void SUPERGREEDY::run(const Graph &G, Stats &stats, int source, int sink, double p, int T)
{
    const auto adj_ = G.get_bidirectional_list();
    int n = G.get_num_nodes();

    std::vector<long long> out_(n, 0); // changed to long long
    std::vector<long long> in_(n, 0);  // changed to long long

    // Initialize out and in degrees for all base nodes.
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

    std::vector<long long> loads(n, 0); // changed to long long
    long long best_cut = -in_[sink];    // changed to long long
    std::cout << "Initial cut: " << best_cut << std::endl;

    for (int t = 0; t < T; t++)
    {
        stats.start_timer();
        std::vector<long long> out = out_;
        std::vector<long long> in = in_;

        std::vector<bool> deleted(n, false);
        deleted[source] = true;
        deleted[sink] = true;
        std::vector<long long> scores(n, 0); // changed to long long

        // Calculate initial marginals as in - out
        for (int i = 0; i < n; ++i)
        {
            scores[i] = loads[i] + (in[i] - out[i]);
        }

        // Add them all to a minheap
        FibHeap heap;
        std::vector<Handle> handles(n);
        for (int i = 0; i < n; ++i)
        {
            handles[i] = heap.push({i, scores[i]});
        }

        long long current_cut = -in[sink]; // changed to long long

        // Start greedy peeling
        while (!heap.empty())
        {
            auto max = heap.top();
            heap.pop();
            const auto i = max.key;
            const auto Di = max.priority;
            if (deleted[i])
            {
                continue;
            }

            const long long delta = in[i] - out[i]; // changed to long long
            loads[i] = Di;

            current_cut -= delta;
            best_cut = std::max(best_cut, current_cut);
            deleted[i] = true;

            // Update scores of neighbors
            for (const auto &[nbi, w] : adj_[i])
            {
                if (deleted[nbi])
                {
                    continue;
                }

                scores[nbi] -= 2LL * w; // ensure 2LL is long long
                out[nbi] += w;
                in[nbi] -= w;
                heap.decrease(handles[nbi], {nbi, scores[nbi]});
            }
        }

        stats.pause_timer();

        auto load_vector_norm = std::transform_reduce(
                                    loads.begin(), loads.end(), 0.0,
                                    std::plus<double>(), [t](long long x) { return std::pow(static_cast<double>(x), 2); }) /
                                ((t + 1) * (t + 1));

        log("Iteration ", t, ": Best cut: ", best_cut);
        stats.push(t, best_cut, load_vector_norm);
        log("Load vector norm: ", load_vector_norm);
    }
    log("Final cut:", best_cut);
}
