#pragma once
#include <chrono>
#include <cstdint>
#include <string>
#include <vector>

struct StatEntry
{
    int iter;
    double best_density;
    double load_norm;
    int64_t elapsed_time_ns;

    StatEntry(const int iter, const double best_density, const double load_norm, const int64_t elapsed_time_ns)
        : iter(iter), best_density(best_density), load_norm(load_norm), elapsed_time_ns(elapsed_time_ns)
    {
    }
    StatEntry() = default;
};

class Stats
{
public:
    /**
     * @brief Construct a new Stats object
     * @param max_entries the maximum number of entries to reserve (set to the number of iterations).
     */
    explicit Stats(const int max_entries);

    /**
     * @brief Start the timer.
     */
    void start_timer();

    /**
     * @brief Pause the timer.
     */
    void pause_timer();

    /**
     * @brief Push a new entry to the stats.
     */
    void push(const int iter, const double best_density, const double load_norm);

    /**
     * @brief Write the stats to a file.
     */
    void dump(const std::string &filename) const;

    const std::vector<double> &norms() const;

    const std::vector<double> &densities() const;

private:
    std::vector<StatEntry> entries_;
    std::vector<double> norms_;
    std::vector<double> densities_;
    int size_;
    bool running_;
    std::chrono::high_resolution_clock::time_point start_time_;
    int64_t elapsed_ns_ = 0;
};