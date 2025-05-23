// pushrelabel.h
#pragma once

#include "graph.h"
#include "stats.h"
#ifdef __cplusplus
extern "C" {
#endif
int run_pushrelabel(int argc, char **argv, Stats &stats, Graph &graph);
#ifdef __cplusplus
}
#endif
