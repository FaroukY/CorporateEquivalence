#include "stats.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

Stats::Stats(int max_entries) : entries_(), size_(0), running_(false), elapsed_ns_(0)
{
    entries_.reserve(max_entries);
    norms_.reserve(max_entries);
    densities_.reserve(max_entries);
}

void Stats::start_timer()
{
    if (!running_)
    {
        start_time_ = std::chrono::high_resolution_clock::now();
        running_ = true;
    }
}

void Stats::pause_timer()
{
    if (running_)
    {
        auto end_time = std::chrono::high_resolution_clock::now();
        elapsed_ns_ += std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time_).count();
        running_ = false;
    }
}

void Stats::push(const int iter, const double best_density, const double load_norm)
{
    if (size_ >= entries_.capacity())
        throw std::runtime_error("Attempt to push too many entries");
    entries_.emplace_back(iter,
                          best_density,
                          load_norm,
                          elapsed_ns_);
    norms_.push_back(load_norm);
    densities_.push_back(best_density);
    ++size_;
    elapsed_ns_ = 0; // Reset elapsed time
}

void Stats::dump(const std::string &filename) const
{
    namespace fs = std::filesystem;

    // Create directories if they do not exist
    fs::path filepath(filename);
    fs::path dir = filepath.parent_path();
    if (!dir.empty() && !fs::exists(dir))
    {
        if (!fs::create_directories(dir))
        {
            throw std::runtime_error("Failed to create directories.");
        }
    }

    std::ofstream file(filename);
    if (!file.is_open())
    {
        throw std::runtime_error("Failed to open file.");
    }

    for (const auto &entry : entries_)
    {
        file << entry.iter << " "
             << entry.best_density << " "
             << entry.load_norm << " "
             << entry.elapsed_time_ns << "\n";
    }
}

const std::vector<double> &Stats::norms() const
{
    return norms_;
}

const std::vector<double> &Stats::densities() const
{
    return densities_;
}