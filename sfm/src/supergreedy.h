#pragma once

#include "graph.h" //Include the graph class
#include "stats.h" //Include the stats class

class SUPERGREEDY
{
public:
    /**
     * @brief Perform the SUPERGREEDY algorithm on the loaded graph
     */
    void run(const Graph &G, Stats &stats, int source, int sink, double p = 1.0, int T = 2);
};