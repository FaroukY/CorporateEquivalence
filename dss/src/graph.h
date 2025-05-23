#pragma once
#include <vector>

class Graph
{
private:
    std::vector<std::vector<std::pair<int, int>>> adjacency_list;
    std::vector<std::vector<int>> bidirectional_list;
    int num_nodes;
    int num_edges;
    int index = 0;

public:
    // Constructor
    Graph(int n, int m);

    // Add an edge (undirected)
    void add_edge(int u, int v);

    // Get the adjacency list
    const std::vector<std::vector<std::pair<int, int>>> &get_adjacency_list() const;
    const std::vector<std::vector<int>> &get_bidirectional_list() const;

    // Get number of nodes
    int get_num_nodes() const;

    // Get number of edges
    int get_num_edges() const;
};
