#include "frankwolfe.h"
#include <algorithm>
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

// Linear Optimization Oracle: given current w (iterate), weights, and neighbors info,
// compute the Frank–Wolfe direction q (and update best_score).
static std::vector<double> linearOracle(
    const std::vector<double> &w,
    const std::vector<double> &weights,
    const std::vector<std::vector<int>> &neighbors,
    int numLeft,
    int numRight,
    double &best_score)
{
    // 1. Sort left‐node indices by the ratio w[i] (descending in cover‐value)
    std::vector<int> order(numLeft);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](int a, int b) {
                  // compare marginal‐gain ratios: w[a]/weights[a] < w[b]/weights[b]
                  return w[a] < w[b];
              });

    // 2. Greedily “cover” right‐nodes once, tracking which are deleted
    std::vector<bool> deleted(numRight, false);
    std::vector<double> q(numLeft, 0.0);
    double remaining = std::accumulate(weights.begin(), weights.end(), 0.0);
    // update best_score using the fully‐covered initial weight
    best_score = std::max(best_score, remaining / numLeft);

    for (int i = 0; i < numLeft; ++i)
    {
        int u = order[i];
        double gain = 0.0;
        // neighbors are stored in a single array of size numLeft+numRight:
        // assume right‐nodes are indexed [numLeft .. numLeft+numRight-1]
        for (int v : neighbors[numRight + u])
        {
            if (!deleted[v])
            {
                gain += weights[v];
                deleted[v] = true;
            }
        }
        q[u] = gain;
        remaining -= gain;
        if (i < numLeft - 1)
            best_score = std::max(best_score, remaining / (numLeft - i - 1));
    }
    return q;
}

void runFW(const std::string &filename, int T, Stats &stats)
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

        // --- LO step: get new direction q and update best_score inside
        auto q = linearOracle(w, weights, neighbors, numLeft, numRight, best_score);

        double gamma = 1.0 / (t + 1);
        for (int i = 0; i < numLeft; ++i)
        {
            w[i] = (1 - gamma) * w[i] + gamma * q[i];
        }

        stats.pause_timer(); // calculate statistics for experiments.
        double loadnorm = 0.0;
        for (int i = 0; i < numLeft; ++i)
            loadnorm += w[i] * w[i];
        loadnorm = std::sqrt(loadnorm);

        stats.push(t, best_score, loadnorm);
        std::cout << "Iteration " << t + 1 << ": " << best_score << std::endl;
    }
}