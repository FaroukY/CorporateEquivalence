#include "graph.h"
#include <cassert>

// Constructor
Graph::Graph(int n, int m) : num_nodes(n), num_edges(m)
{
    adjacency_list.resize(n);
    bidirectional_list.resize(n);
}

// Add an edge to the graph
void Graph::add_edge(int u, int v)
{
    assert(u != v);
    adjacency_list[u].push_back({v, index + 1});
    adjacency_list[v].push_back({u, index});
    bidirectional_list[u].push_back(v);
    bidirectional_list[v].push_back(u);
    index += 2;
}

const std::vector<std::vector<std::pair<int, int>>> &Graph::get_adjacency_list() const
{
    return adjacency_list;
}

const std::vector<std::vector<int>> &Graph::get_bidirectional_list() const
{
    return bidirectional_list;
}

// Get number of nodes
int Graph::get_num_nodes() const
{
    return num_nodes;
}

// Get number of edges
int Graph::get_num_edges() const
{
    return num_edges;
}