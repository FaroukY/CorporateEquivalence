// flow.h
#pragma once

#include "stats.h"
#include <string>

#ifdef __cplusplus
extern "C" {
#endif
int runFlow(std::string &filename, Stats &stats);
#ifdef __cplusplus
}
#endif
