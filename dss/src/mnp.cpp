#include "mnp.h"
#include "graph.h"
#include "logger.h"
#include "stats.h"
#include <execution>
#include <numeric>
#include <tbb/tbb.h>

namespace {
    double get_numerator(const std::vector<double> &y, const std::vector<double> &degrees, double p)
    {
        const auto degree_sum = std::transform_reduce(
            std::execution::par, degrees.begin(), degrees.end(), 0.0,
            std::plus<double>(), [p](double x) { return std::pow(x, p); });
        const auto y_sum = std::accumulate(
            y.begin(), y.end(), 0.0);
        return degree_sum - 2 * y_sum;
    }

    Eigen::VectorXd LO(const std::vector<double> &y, double p, double numerator, std::vector<double> degrees, const std::vector<std::vector<int>> &adj,
                       const Eigen::VectorXd &loads)
    {
        int n = static_cast<int>(adj.size());
        Eigen::VectorXd q = Eigen::VectorXd::Zero(n);
        std::vector<int> indices(n);
        std::iota(indices.begin(), indices.end(), 0);

        // Sort node indices in ascending order of their current load
        std::sort(std::execution::par, indices.begin(), indices.end(), [&loads](int a, int b) { return loads(a) < loads(b); });

        std::vector<bool> deleted(n, false);

        // Set q(vi) = f(Si) - f(Si-1)
        for (int i = 0; i < indices.size(); i++)
        {
            auto idx = indices[i];
            // log("Peeling node ", idx, " with degree ", degrees[idx], " and numerator ", numerator);
            auto numerator_delta = std::pow(degrees[idx], p) - 2 * y[idx];
            for (auto nbi = 0; nbi < adj[idx].size(); nbi++)
            {
                const auto nb = adj[idx][nbi];
                if (deleted[nb])
                    continue; // Skip deleted nodes

                if (degrees[nb] > 0)
                    numerator_delta += std::pow(degrees[nb], p) - std::pow(degrees[nb] - 1, p);
                degrees[nb]--;
            }
            q(idx) = numerator_delta;

            deleted[idx] = true;
            degrees[idx] = 0;
        }
        return q;
    }

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
} // namespace

// Core Min-Norm-Point algorithm implementation
// Returns final point x, its squared norm, and history of norms per iteration
std::tuple<Eigen::VectorXd, double>
MNP::min_norm_alg(Stats &stats, const std::vector<double> &y, double p, double numerator, const std::vector<double> &degrees, const std::vector<std::vector<int>> &adj,
                  int n, int m, int iterations, bool is_contra)
{

    // 1) Initialize x = LO(adj, zero_vector)
    Eigen::VectorXd x = LO(y, p, numerator, degrees, adj, Eigen::VectorXd::Zero(n));
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
        Eigen::VectorXd q = LO(y, p, numerator, degrees, adj, x);
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
        best_density_ = std::max(best_density_, extract_density(indices, y, p, numerator, degrees, adj, n, m, is_contra));

        log("Iteration ", t, ": ");
        log("Norm value: ", load_vector_norm);
        log("Best density: ", best_density_);
        log("Max b: ", best_x.maxCoeff(), ", Min b: ", best_x.minCoeff());

        stats.push(t, best_density_, load_vector_norm);
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
    const auto best_density = extract_density(indices, y, p, numerator, degrees, adj, n, m, is_contra);
    // 8) Print minimum norm over iterations
    double min_norm_val = *std::min_element(std::execution::par, stats.norms().begin(), stats.norms().end());

    log("Final stats: ");
    log("Min Norm value: ", min_norm_val);
    log("Final density: ", best_density);

    return {x, x.squaredNorm()};
}

void MNP::run(const Graph &G, Stats &stats, const std::vector<double> &y, double p, int T, bool is_contra)
{
    const auto adj_ = G.get_bidirectional_list();
    int n = G.get_num_nodes();
    int m = G.get_num_edges();
    log("Running MNP on graph with ", n, " nodes and ", m, " edges.");

    std::vector<double> degrees_(n, 0);
#pragma omp parallel for
    for (int i = 0; i < n; ++i)
    {
        degrees_[i] = adj_.at(i).size();
    }
    const auto numerator_ = get_numerator(y, degrees_, p);

    // Run the min-norm algorithm for T iterations and capture outputs
    auto [x, norm] = min_norm_alg(stats, y, p, numerator_, degrees_, adj_, n, m, T, is_contra);
    log("Final norm squared: ", norm);
}

double MNP::extract_density(const std::vector<int> &indices,
                            const std::vector<double> &y,
                            double p,
                            double numerator,
                            std::vector<double> degrees,
                            const std::vector<std::vector<int>> &adj,
                            const int n, const int m, bool is_contra) const
{
    double highest_density = std::pow(numerator / n, 1.0 / p);
    std::vector<bool> deleted(n, false);

    auto denominator = n;
    auto denom = is_contra ? 1 : n;
    for (int i = 0; i < indices.size(); i++)
    {
        auto idx = indices[i];
        auto numerator_delta = -std::pow(degrees[idx], p) + 2 * y[idx];
        for (auto nbi = 0; nbi < adj[idx].size(); nbi++)
        {
            const auto nb = adj[idx][nbi];
            if (deleted[nb])
                continue; // Skip deleted nodes

            if (degrees[nb] > 0)
                numerator_delta -= std::pow(degrees[nb], p) - std::pow(degrees[nb] - 1, p);
            degrees[nb]--;
        }
        numerator = numerator + numerator_delta;
        denominator--;
        deleted[idx] = true;
        degrees[idx] = 0;

        denom = is_contra ? 1 : std::max(denominator, 1);
        highest_density = std::max(highest_density,
                                   std::pow(numerator / denom, 1.0 / p));
    }

    return highest_density;
}