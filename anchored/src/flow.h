// flow.h
#pragma once

#include "graph.h"
#include "stats.h"

#ifdef __cplusplus
extern "C" {
#endif
int run_flow(const Graph &G, Stats &stats, const std::vector<int> &penalty);
#ifdef __cplusplus
}
#endif
