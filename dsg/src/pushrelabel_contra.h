// pushrelabel.h
#pragma once

#include "stats.h"

#ifdef __cplusplus
extern "C" {
#endif
int run_pushrelabel_contra(int argc, char **argv, Stats &stats, const std::string &outfile_name);
#ifdef __cplusplus
}
#endif
