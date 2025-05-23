#include <algorithm>
#include <cmath>
#include <execution>
#include <iostream>
#include <limits>
#include <numeric>
#include <tbb/tbb.h>
#include <vector>

#include "logger.h"
#include "mnp.h"

namespace {
    // LO: Linear optimization oracle for the Min-Norm Point (MNP) algorithm
    // This is essentially Edmond's algorithm (the Greedy Algorithm)
    // Given the current load vector, computes q = LO(G, n, loads)
    // by sorting nodes by load in ascending order v1, v2, ..., vn
    // then letting Si = {v1, v2, ..., vi} and q(vi) = f(Si)-f(Si-1) and return q.
    // This corresponds to deleting the nodes in order v1, ..., vn, then for each i,
    // set q[vi] as the number of edges left connected to vi.
    Eigen::VectorXd LO(const std::vector<std::vector<int>> &adj,
                       const Eigen::VectorXd &loads)
    {
        int n = static_cast<int>(adj.size());
        Eigen::VectorXd q = Eigen::VectorXd::Zero(n);
        std::vector<int> indices(n);
        std::iota(indices.begin(), indices.end(), 0);

        // Sort node indices in ascending order of their current load
        std::sort(indices.begin(), indices.end(), [&loads](int a, int b) { return loads(a) < loads(b); });

        // Track which nodes are "deleted" (processed) to count only active edges
        std::vector<bool> deleted(n, false);
        for (int idx : indices)
        {
            // For each neighbor of idx not yet deleted, increment q(idx)
            for (int neighbor : adj[idx])
            {
                if (!deleted[neighbor])
                    q(idx)++;
            }
            // Mark this node as deleted so later nodes don't count its edges
            deleted[idx] = true;
        }
        return q;
    }

    // AffineMinimizer: finds the affine combination of points in S that minimizes norm
    // Uses an LLT (Cholesky) solve to avoid explicit matrix inversion
    std::pair<Eigen::VectorXd, Eigen::VectorXd>
    AffineMinimizer(const std::vector<Eigen::VectorXd> &S)
    {
        int n = S[0].size();                // dimension of each vector
        int k = static_cast<int>(S.size()); // number of extreme points
        Eigen::MatrixXd B(n, k);
        for (int j = 0; j < k; j++)
        {
            B.col(j) = S[j]; // build matrix B whose columns are the points
        }
        // M = Bᵀ * B
        Eigen::MatrixXd M = B.transpose() * B;
        Eigen::VectorXd ones = Eigen::VectorXd::Ones(k);
        // Solve M * sol = 1 via Cholesky factorization
        Eigen::LLT<Eigen::MatrixXd> llt(M);
        Eigen::VectorXd sol = llt.solve(ones);
        // Normalize to enforce affine constraint (coefficients sum to 1)
        double denominator = sol.sum();
        Eigen::VectorXd alpha = sol / denominator;
        // Compute the minimizer y = B * alpha
        Eigen::VectorXd y = B * alpha;
        return {y, alpha};
    }

} // end anonymous namespace

// Core Min-Norm-Point algorithm implementation
// Returns final point x, its squared norm, and history of norms per iteration
std::tuple<Eigen::VectorXd, double>
MNP::min_norm_alg(Stats &stats, const std::vector<std::vector<int>> &adj,
                  int n, int m, int iterations)
{
    // 1) Initialize x = LO(adj, zero_vector)
    Eigen::VectorXd x = LO(adj, Eigen::VectorXd::Zero(n));
    Eigen::VectorXd best_x = x;
    double best_heur = std::numeric_limits<double>::infinity();

    // Maintain set S of extreme points and corresponding lambdas
    std::vector<Eigen::VectorXd> S{x};
    std::vector<double> lambdas{1.0};

    // Sorted indices for tracking best density
    std::vector<int> indices(n);
    std::iota(indices.begin(), indices.end(), 0);

    // Main loop: perform 'iterations' rounds of MNP
    for (int t = 0; t < iterations; t++)
    {
        stats.start_timer();
        // 2) Compute new extreme direction q = LO(adj, x)
        Eigen::VectorXd q = LO(adj, x);
        // 3) Heuristic: xᵀx - xᵀq to track progress
        double current_heur = x.squaredNorm() - x.dot(q);
        if (current_heur < best_heur)
        {
            best_heur = current_heur;
            best_x = x;
        }

        // 4) Add q to the set S, lambda for new point initialized to 0
        S.push_back(q);
        lambdas.push_back(0.0);

        // 5) Affine minimization with nonnegativity enforcement
        while (true)
        {
            // Compute minimizer y and coefficients alpha
            auto [y, alpha] = AffineMinimizer(S);
            // If all alpha are nonnegative, accept y as new x
            if (alpha.minCoeff() >= 0.0)
            {
                x = y;
                lambdas.assign(alpha.data(), alpha.data() + alpha.size());
                break;
            }
            // Otherwise, perform line search: find smallest theta that makes lambdas+theta*alpha >= 0
            double theta = std::numeric_limits<double>::infinity();
            for (int j = 0; j < alpha.size(); ++j)
            {
                if (alpha(j) < 0.0)
                {
                    double candidate = lambdas[j] / (lambdas[j] - alpha(j));
                    theta = std::min(theta, candidate);
                }
            }
            // Update x and lambdas via convex combination
            x = theta * y + (1.0 - theta) * x;
            for (int j = 0; j < lambdas.size(); ++j)
            {
                lambdas[j] = (1.0 - theta) * lambdas[j] + theta * alpha(j);
            }
            // Remove any points with nonpositive lambda (inactive extremes)
            for (int j = lambdas.size() - 1; j >= 0; --j)
            {
                if (lambdas[j] <= 0.0)
                {
                    lambdas.erase(lambdas.begin() + j);
                    S.erase(S.begin() + j);
                }
            }
        }
        // Calculate and store stats
        stats.pause_timer();
        std::sort(std::execution::par, indices.begin(), indices.end(), [&best_x](int a, int b) { return best_x(a) < best_x(b); });
        const auto load_vector_norm = x.squaredNorm();
        const auto best_density = extract_density(indices, adj, n, m);

        log("Iteration ", t, ": ");
        log("Norm value: ", load_vector_norm);
        log("Best density: ", best_density);
        log("Max b: ", best_x.maxCoeff(), ", Min b: ", best_x.minCoeff());

        stats.push(t, best_density, load_vector_norm);
    }

    // 6) Debug: print first few coords of best_x and norm history
    for (int i = 0; i < std::min(10, n); ++i)
        log("Best x[", i, "] = ", best_x(i));
    log("");
    for (const double norm : stats.norms())
        log(norm, ", ");
    log("");
    log("Best x sum: ", best_x.sum());

    // 7) Fractional peeling on best_x to compute final density
    const auto best_density = extract_density(indices, adj, n, m);
    // 8) Print minimum norm over iterations
    double min_norm_val = *std::min_element(std::execution::par, stats.norms().begin(), stats.norms().end());

    log("Final stats: ");
    log("Min Norm value: ", min_norm_val);
    log("Final density: ", best_density);

    return {x, x.squaredNorm()};
}

// Ties the MNP machinery into the Graph interface.
void MNP::run(const Graph &graph, Stats &stats, int T)
{
    int n = graph.get_num_nodes();
    int m = graph.get_num_edges();
    log("Running MNP on graph with ", n, " nodes and ", m, " edges.");
    // Use undirected neighbor list for density calculation
    std::vector<std::vector<int>> adjacency_list = graph.get_bidirectional_list();

    // Run the min-norm algorithm for T iterations and capture outputs
    auto [x, norm] = min_norm_alg(stats, adjacency_list, n, m, T);
    log("Final norm squared: ", norm);
}

double MNP::extract_density(const std::vector<int> &indices,
                            const std::vector<std::vector<int>> &adj,
                            const int n, const int m) const
{
    double highest_density = static_cast<double>(m) / n;
    int edges_left = m;
    int num_nodes = n;
    std::vector<bool> visited(n, false);
    for (int idx : indices)
    {
        visited[idx] = true;
        for (int nei : adj[idx])
            if (!visited[nei])
                edges_left--;
        num_nodes--;
        if (num_nodes > 0)
            highest_density = std::max(highest_density,
                                       static_cast<double>(edges_left) / num_nodes);
    }
    return highest_density;
}