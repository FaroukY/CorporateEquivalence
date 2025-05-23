#include "graph.h"
#include <cassert>

// Constructor
Graph::Graph(int n, int m) : num_nodes(n), num_edges(m)
{
    adjacency_list.resize(n);
    bidirectional_list.resize(n);
    edgeMap.resize(2 * m);
}

// Add an edge to the graph
void Graph::add_edge(int u, int v)
{
    assert(u != v);
    adjacency_list[v].push_back({u, index});
    edgeMap[index] = v;
    adjacency_list[u].push_back({v, index + 1});
    edgeMap[index + 1] = u;
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

const std::vector<int> &Graph::get_edge_map() const
{
    return edgeMap;
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