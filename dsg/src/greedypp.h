#pragma once

#include "graph.h"
#include "stats.h"

class GreedyPeelingPP
{
public:
    void run(const Graph &graph, Stats &stats, int iterations = 1);
    double get_density() const
    {
        return mm_density;
    }
    std::vector<int> get_ans() const
    {
        return w;
    }

private:
    double mm_density = 0;
    std::vector<int> w;
};
