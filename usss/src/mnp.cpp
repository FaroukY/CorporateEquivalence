#include <algorithm>
#include <cmath>
#include <execution>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include "mnp.h"
#include "stats.h"

namespace {

    // Linear‐Optimization Oracle for USSS
    Eigen::VectorXd LO(
        const Eigen::VectorXd &w,
        const std::vector<double> &weights,
        const std::vector<std::vector<int>> &neighbors,
        int numLeft,
        int numRight,
        double &best_score)
    {
        Eigen::VectorXd q = Eigen::VectorXd::Zero(numLeft);

        // sort left‐nodes by current load w
        std::vector<int> order(numLeft);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(),
                  [&](int a, int b) { return w[a] < w[b]; });

        // greedily cover right‐nodes
        std::vector<bool> deleted(numRight, false);
        double remaining = std::accumulate(weights.begin(), weights.end(), 0.0);

        for (int i = 0; i < numLeft; ++i)
        {
            int u = order[i];
            double gain = 0.0;
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

    // Affine minimizer via Cholesky, exactly as in DSG‐MNP.
    std::pair<Eigen::VectorXd, Eigen::VectorXd>
    AffineMinimizer(const std::vector<Eigen::VectorXd> &S)
    {
        int n = S[0].size();
        int k = static_cast<int>(S.size());
        Eigen::MatrixXd B(n, k);
        for (int j = 0; j < k; ++j)
        {
            B.col(j) = S[j];
        }
        Eigen::MatrixXd M = B.transpose() * B;
        Eigen::LLT<Eigen::MatrixXd> llt(M);
        Eigen::VectorXd sol = llt.solve(Eigen::VectorXd::Ones(k));
        Eigen::VectorXd alpha = sol / sol.sum();
        return {B * alpha, alpha};
    }

    // Core MNP loop for USSS
    void min_norm_alg_usss(
        Stats &stats,
        const std::vector<std::vector<int>> &neighbors,
        const std::vector<double> &weights,
        int numLeft,
        int numRight,
        int iterations)
    {
        double best_score = 0;

        // 1) initialize x = LO_usss(0)
        Eigen::VectorXd x = Eigen::VectorXd::Zero(numLeft);
        x = LO(x, weights, neighbors, numLeft, numRight, best_score);
        Eigen::VectorXd best_x = x;
        double best_heur = std::numeric_limits<double>::infinity();

        // maintain extreme set S and lambdas
        std::vector<Eigen::VectorXd> S{x};
        std::vector<double> lambdas{1.0};

        for (int t = 0; t < iterations; ++t)
        {
            stats.start_timer();

            // 2) oracle
            Eigen::VectorXd q = LO(x, weights, neighbors, numLeft, numRight, best_score);

            // 3) heuristic tracking
            double curr = x.squaredNorm() - x.dot(q);
            if (curr < best_heur)
            {
                best_heur = curr;
                best_x = x;
            }

            // 4) add new extreme
            S.push_back(q);
            lambdas.push_back(0.0);

            // 5) projection loop
            while (true)
            {
                auto [y, alpha] = AffineMinimizer(S);
                if (alpha.minCoeff() >= 0.0)
                {
                    x = y;
                    lambdas.assign(alpha.data(), alpha.data() + alpha.size());
                    break;
                }
                double theta = std::numeric_limits<double>::infinity();
                for (int j = 0; j < alpha.size(); ++j)
                {
                    if (alpha[j] < 0.0)
                    {
                        double cand = lambdas[j] / (lambdas[j] - alpha[j]);
                        theta = std::min(theta, cand);
                    }
                }
                x = theta * y + (1 - theta) * x;
                for (int j = 0; j < (int)lambdas.size(); ++j)
                {
                    lambdas[j] = (1 - theta) * lambdas[j] + theta * alpha[j];
                }
                // drop inactive extremes
                for (int j = (int)lambdas.size() - 1; j >= 0; --j)
                {
                    if (lambdas[j] <= 0.0)
                    {
                        lambdas.erase(lambdas.begin() + j);
                        S.erase(S.begin() + j);
                    }
                }
            }
            stats.pause_timer();
            double normVal = x.squaredNorm();
            normVal = sqrt(normVal);
            stats.push(t, best_score, normVal);
        }
    }

} // anonymous namespace

void runMNP(const std::string &filename, int T, Stats &stats)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error opening file " << filename << "\n";
        return;
    }

    int numLeft, numRight;
    file >> numLeft >> numRight;

    std::vector<double> weights(numRight);
    for (int i = 0; i < numRight; ++i)
    {
        file >> weights[i];
    }

    std::vector<std::vector<int>> neighbors(numLeft + numRight);
    int u, v;
    while (file >> u >> v)
    {
        neighbors[v].push_back(u);
        neighbors[u].push_back(v);
    }
    file.close();
    min_norm_alg_usss(stats, neighbors, weights, numLeft, numRight, T);
}
