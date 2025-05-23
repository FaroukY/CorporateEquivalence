#include "supergreedy.h"
#include <algorithm>
#include <boost/heap/fibonacci_heap.hpp>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

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

void runIP(const std::string &filename, int T, Stats &stats)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error opening file " << filename << std::endl;
        return;
    }

    int numLeft, numRight;
    file >> numLeft >> numRight;

    std::vector<double> weights(numRight);
    for (int i = 0; i < numRight; ++i)
    {
        file >> weights[i];
    }

    std::vector<std::vector<int>> neighbors(numLeft + numRight, std::vector<int>());
    int u, v;
    while (file >> u >> v)
    {
        neighbors[v].push_back(u);
        neighbors[u].push_back(v);
    }
    file.close();

    double best_score = 0.0;
    std::vector<double> w(numLeft, 0.0);

    for (int t = 0; t < T; ++t)
    {
        stats.start_timer();

        // Build degree-based priority queue
        FibHeap pq;
        std::vector<Handle> handles(numLeft);

        std::vector<double> loads(numLeft, 0.0);
        for (int u = 0; u < numLeft; ++u)
        {
            double score = 0.0;
            for (int v : neighbors[numRight + u])
                score += weights[v];
            loads[u] = w[u] + score;
            handles[u] = pq.push({u, loads[u]});
        }

        std::vector<bool> deleted(numRight, false);
        std::unordered_set<int> processed;
        double remaining = std::accumulate(weights.begin(), weights.end(), 0.0);
        best_score = std::max(best_score, remaining / numLeft);
        while (!pq.empty())
        {
            auto top = pq.top();
            pq.pop();
            int u = top.key;

            if (processed.count(u))
                continue;
            processed.insert(u);

            double score = 0;
            for (int v : neighbors[numRight + u])
            {
                if (!deleted[v])
                {
                    score += weights[v];
                    deleted[v] = true;
                    // Update v's neighbors loads
                    for (auto u_neighbor : neighbors[v])
                    {
                        int u_neighbor_idx = u_neighbor - numRight;
                        if (!processed.count(u_neighbor_idx))
                        {
                            loads[u_neighbor_idx] -= weights[v];
                            pq.decrease(handles[u_neighbor_idx], {u_neighbor_idx, loads[u_neighbor_idx]});
                        }
                    }
                }
            }
            remaining -= score;
            if (processed.size() < numLeft)
                best_score = std::max(best_score, remaining / (numLeft - processed.size()));
            w[u] += score; // accumulate count into w[u]
        }
        stats.pause_timer(); // calculate statistics for experiments.

        double loadnorm = 0.0;
        for (int i = 0; i < numLeft; ++i)
            loadnorm += w[i] * w[i] / ((t + 1) * (t + 1));
        loadnorm = std::sqrt(loadnorm);

        stats.push(t, best_score, loadnorm);
        std::cout << "Iteration " << t + 1 << ": " << best_score << std::endl;
    }
}