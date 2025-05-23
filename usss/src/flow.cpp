// High-performance push-relabel based exact DSP solver
extern "C" {
#include "external/exactDSP-cpp/hi_pr.h"
}

// Standard library includes for IO, containers, timing, etc.
#include "flow.h"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
using namespace std;

// INF: large capacity to model infinite edges in the flow network
const int INF = (int)1e9;

// add_arc: add directed arc u->v with given capacity, plus reverse arc v->u of zero capacity
void add_arc(int u, int v, int capacity, node *nodes, arc *arcs, cType *cap, long *cur_arc)
{
    arcs[cur_arc[u]].head = nodes + v;        // set head pointer for forward arc
    arcs[cur_arc[u]].rev = arcs + cur_arc[v]; // link to reverse arc
    arcs[cur_arc[v]].head = nodes + u;        // set head for reverse arc
    arcs[cur_arc[v]].rev = arcs + cur_arc[u]; // link reverse to forward
    cap[cur_arc[u]] = capacity;               // forward capacity
    cap[cur_arc[v]] = 0;                      // reverse capacity = 0
    ++cur_arc[u];                             // advance write pointer for u
    ++cur_arc[v];                             // advance write pointer for v
}

// nontrivial: returns true if more than one node remains on the source side of the cut
bool nontrivial(int n_nodes, node *nodes)
{
    int res = 0;
    for (int i = 0; i < n_nodes; ++i)
    {
        if (nodes[i].d < n_nodes) // node label < n_nodes indicates reachability
            ++res;
    }
    return res > 1;
}

// run_pushrelabel: main entry, constructs the flow network and iteratively extracts densest subgraph
extern "C" int runFlow(std::string &filename, Stats &stats)
{
    // ACCURACY NUMBER DETERMINES HOW MANY DECIMALS WE COMPUTE THINGS TO
    // HI_PR ONLY TAKES INTEGER CAPACITIES, so we scale all capacities by ACCURACY
    auto startio = chrono::steady_clock::now();

    int ACCURACY = 100; // Accuracy 1/100 = 0.01

    // --- 1) Read input ---
    ifstream file(filename);
    if (!file.is_open())
    {
        cerr << "Error opening file " << filename << "\n";
        return -1;
    }

    int numLeft, numRight;
    file >> numLeft >> numRight;

    // real‐valued weights on the right side
    vector<double> w_double(numRight);
    for (int j = 0; j < numRight; ++j)
        file >> w_double[j];

    // Build neighbors exactly as you did (global indexing):
    //  Right: 0..numRight-1    Left: numRight..numRight+numLeft-1
    vector<vector<int>> neighbors(numLeft + numRight);
    int u, v;
    while (file >> u >> v)
    {
        neighbors[u].push_back(v);
        neighbors[v].push_back(u);
    }
    file.close();

    // Convert weights to integer capacities
    vector<int> w_int(numRight);
    for (int j = 0; j < numRight; ++j)
        w_int[j] = int(round(w_double[j] * ACCURACY));

    // Read input sizes
    int n = numLeft;  // number of left vertices
    int m = numRight; // number of right vertices

    auto start = chrono::steady_clock::now();

    // Define special nodes and total node count
    int source = n + m;
    int sink = n + m + 1;
    int n_nodes = n + m + 2;

    // deg array will hold prefix sums of arc counts per node
    long *deg = new long[n_nodes];
    memset(deg, 0, sizeof(long) * n_nodes);

    // Initialize out-degrees for prefix calculation
    deg[source] = n; // source connects to all left vertices
    deg[sink] = m;   // all rightvertices connect to sink

    for (int u = 0; u < numLeft; ++u)
        deg[numRight + u] = 1 + neighbors[numRight + u].size(); // 1 from s, and rest from right neighbors

    for (int v = 0; v < numRight; ++v)
        deg[v] = 1 + neighbors[v].size(); // 1 from sink, and rest from left neighbors

    for (int i = 1; i < n_nodes; ++i)
        deg[i] += deg[i - 1]; // prefix sum

    int n_arcs = deg[n_nodes - 1]; // total number of directed arcs

    // Allocate write pointers and flow network structures
    long *cur_arc = new long[n_nodes];
    node *nodes_ = new node[n_nodes + 1];
    arc *arcs = new arc[n_arcs];
    cType *cap = new cType[n_arcs];

    // Initialize per-node arc start pointers
    for (int i = 0; i < n_nodes; ++i)
    {
        cur_arc[i] = (i == 0 ? 0 : deg[i - 1]);
        nodes_[i].first = arcs + cur_arc[i];
    }

    // Add arcs: source -> left vertices (0 capacity initially)
    for (int u = 0; u < numLeft; ++u)
        add_arc(source, numRight + u, 0, nodes_, arcs, cap, cur_arc);

    // Add arcs: rightNodes -> sink (w_int[i] capacity)
    for (int v = 0; v < numRight; ++v)
        add_arc(v, sink, w_int[v], nodes_, arcs, cap, cur_arc);

    // Add arcs: numLeft -> numRight (INF capacity)
    for (int u = 0; u < numLeft; ++u)
        for (int v : neighbors[numRight + u])
            add_arc(numRight + u, v, INF, nodes_, arcs, cap, cur_arc);

    // Initialize subgraph membership (all vertices included)
    vector<bool> subg(numLeft, true);
    node *j;

    // Initial density: w(V)/|N(V)| = sum_w / n
    double sumW = accumulate(w_int.begin(), w_int.end(), 0.0);
    double lambda = (sumW / n);
    int c = (int)lambda, prev_c = -1;
    std::cout << c << " " << lambda << std::endl;
    // Iteratively refine candidate subgraph until density stabilizes
    int iteration = 0;

    while (c != prev_c)
    {
        stats.start_timer();
        prev_c = c;

        // Update capacities on source->left vertex arcs based on subg membership
        for (int u = 0; u < numLeft; ++u)
            cap[nodes_[source].first - arcs + u] = subg[u] ? c : 0;

        // Update capacities on edge->sink arcs based on both endpoints in subg
        for (int v = 0; v < numRight; ++v)
        {
            bool ok = true;
            for (int u : neighbors[v])
                if (!subg[u - numRight])
                    ok = false;
            long start_c = nodes_[v].first - arcs;
            cap[start_c] = ok ? w_int[v] : 0;
        }

        // Run push-relabel min_cut on the current network
        min_cut(n_nodes, n_arcs / 2,
                nodes_, arcs, cap,
                nodes_ + source, nodes_ + sink, 0);

        // Extract new subgraph: vertices still reachable from source
        fill(subg.begin(), subg.end(), false);
        int new_n = 0;

        if (nontrivial(n_nodes, nodes_))
        {
            forAllNodes(j) if (j->d < n_nodes && nNode(j) >= numRight && nNode(j) < numLeft + numRight && cap[nodes_[source].first - arcs + nNode(j) - numRight] > 0)
            {
                subg[nNode(j) - numRight] = true;
                new_n++;
            }
        }

        std::cout << "new_n: " << new_n << std::endl;
        // Recompute density c for the updated subgraph
        double newW = 0.0;
        for (int v = 0; v < numRight; ++v)
        {
            bool ok = true;
            for (int u : neighbors[v])
                if (!subg[u - numRight])
                    ok = false;
            if (ok)
                newW += w_int[v];
        }
        c = static_cast<int>((1.0 * newW / max(1, new_n)));
        const auto current_density = (double)c / ACCURACY;
        stats.pause_timer();
        cerr << "Current density estimate: " << current_density << endl;

        stats.push(iteration++, current_density, 0.0 /* load_norm not applicable */);
    }
    delete[] deg;
    delete[] cur_arc;
    delete[] nodes_;
    delete[] arcs;
    delete[] cap;

    return 0;
}