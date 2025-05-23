#pragma once
#include "graph.h"
#include <cstdio>  // Needed for `fread()`
#include <cstdlib> // Needed for `EOF`
#include <string>

class Parser
{
public:
    Graph read_graph_from_structured_file(const std::string &filename);
};
