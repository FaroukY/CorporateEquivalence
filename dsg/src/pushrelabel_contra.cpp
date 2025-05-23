// High-performance push-relabel based exact DSP solver
extern "C" {
#include "../../external/exactDSP-cpp/hi_pr.h"
}

// Standard library includes for IO, containers, timing, etc.
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

#include "pushrelabel_contra.h"
using namespace std;

namespace {
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

    // GET_CHAR: buffered character reader for fast input
    inline char GET_CHAR()
    {
        const int maxn = 131072;
        static char buf[maxn], *p1 = buf, *p2 = buf;
        return p1 == p2 && (p2 = (p1 = buf) + fread(buf, 1, maxn, stdin), p1 == p2)
                   ? EOF
                   : *p1++;
    }

    // getInt: parse next integer from stdin using GET_CHAR
    inline int getInt()
    {
        int res = 0;
        char c = GET_CHAR();
        while (c < '0')
            c = GET_CHAR();
        while (c >= '0')
        {
            res = res * 10 + (c - '0');
            c = GET_CHAR();
        }
        return res;
    }
} // namespace

// run_pushrelabel: main entry, constructs the flow network and iteratively extracts densest subgraph
extern "C" int run_pushrelabel_contra(int argc, char **argv, Stats &stats, const std::string &outfile_name)
{
    // ACCURACY NUMBER DETERMINES HOW MANY DECIMALS WE COMPUTE THINGS TO
    // HI_PR ONLY TAKES INTEGER CAPACITIES, so we scale all capacities by ACCURACY
    auto startio = chrono::steady_clock::now();

    int ACCURACY = atoi(argv[1]); // scale factor for densities

    // Read input sizes
    int n = getInt(); // number of vertices
    int m = getInt(); // number of edges
    int k = 2;        // edge size (edges case)

    // Read all edge endpoints
    int *edges = new int[m * k];
    for (int i = 0; i < m * k; ++i)
        edges[i] = getInt();

    auto endio = chrono::steady_clock::now();
    cout << "Time for reading input: "
         << chrono::duration_cast<chrono::milliseconds>(endio - startio).count()
         << " ms" << endl;

    auto start = chrono::steady_clock::now();

    // Define special nodes and total node count
    int source = n + m;
    int sink = n + m + 1;
    int n_nodes = n + m + 2;

    // deg array will hold prefix sums of arc counts per node
    long *deg = new long[n_nodes];

    // Initialize out-degrees for prefix calculation
    deg[source] = n; // source connects to all vertices
    deg[sink] = m;   // all edges connect to sink
    for (int i = 0; i < n; ++i)
        deg[i] = 1;
    for (int i = 0; i < m * k; ++i)
        ++deg[edges[i]];
    for (int i = 0; i < m; ++i)
        deg[n + i] = k + 1;
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

    // Add arcs: source -> vertices (0 capacity)
    for (int i = 0; i < n; ++i)
        add_arc(source, i, 0, nodes_, arcs, cap, cur_arc);

    // Add arcs: edge-nodes -> sink (ACCURACY capacity)
    for (int i = 0; i < m; ++i)
        add_arc(n + i, sink, ACCURACY, nodes_, arcs, cap, cur_arc);

    // Add arcs: vertices -> their incident edge-nodes (infinite capacity)
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < k; ++j)
            add_arc(edges[k * i + j], n + i, INF, nodes_, arcs, cap, cur_arc);

    // Initialize subgraph membership (all vertices included)
    vector<bool> subg(n, true);
    node *j;

    int prev_c = -1;
    int c = static_cast<int>((1.0 * m / n) * ACCURACY); // initial density estimate

    // Iteratively refine candidate subgraph until density stabilizes
    int iteration = 0;
    while (c != prev_c)
    {
        stats.start_timer();
        prev_c = c;

        // Update capacities on source->vertex arcs based on subg membership
        for (int u = 0; u < n; ++u)
            cap[nodes_[source].first - arcs + u] = subg[u] ? c : 0;

        // Update capacities on edge->sink arcs based on both endpoints in subg
        for (int e = 0; e < m; ++e)
        {
            long start_c = nodes_[n + e].first - arcs;
            bool keep = subg[edges[2 * e]] && subg[edges[2 * e + 1]];
            cap[start_c] = keep ? ACCURACY : 0;
        }

        // Run push-relabel min_cut on the current network
        min_cut(n_nodes, n_arcs / 2,
                nodes_, arcs, cap,
                nodes_ + source, nodes_ + sink, 0);

        // Extract new subgraph: vertices still reachable from source
        fill(subg.begin(), subg.end(), false);
        if (nontrivial(n_nodes, nodes_))
        {
            forAllNodes(j) if (j->d < n_nodes && nNode(j) < n && cap[nodes_[source].first - arcs + nNode(j)] > 0)
                subg[nNode(j)] = true;
        }

        // Recompute density c for the updated subgraph
        stats.pause_timer();
        int subg_size = accumulate(subg.begin(), subg.end(), 0);
        int subg_edges = 0;
        for (int i = 0; i < m; ++i)
            subg_edges += (subg[edges[i * k]] && subg[edges[i * k + 1]]);
        c = static_cast<int>((1.0 * subg_edges / max(1, subg_size)) * ACCURACY);
        const auto current_density = (double)c / ACCURACY;
        cerr << "Current density estimate: " << current_density << " " << subg_edges << " " << subg_size
             << endl;

        stats.push(iteration++, current_density, 0.0 /* load_norm not applicable */);
    }

    // Save the final subgraph
    std::cout << "Saving final subgraph to " << outfile_name << std::endl;
    std::ofstream outfile(outfile_name, std::ios::app); // append mode
    for (int i = 0; i < n; ++i)
    {
        if (subg[i])
            outfile << i << " ";
    }
    outfile << std::endl;
    // Alongside density
    outfile << (double)c / ACCURACY << std::endl;
    outfile.close();

    // Output total time for network solve and final density
    auto end = chrono::steady_clock::now();
    cout << "Time for finding solution: "
         << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << " ms" << endl;
    cerr << "Density is " << (double)c / ACCURACY << endl;

    return 0;
}
