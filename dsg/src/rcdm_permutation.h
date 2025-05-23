#pragma once

#include "graph.h"
#include "stats.h"

void runRCDM(const Graph &G, int T, Stats &stats, bool provided_best_loads = false);
