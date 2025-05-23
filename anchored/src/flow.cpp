// flow.cpp

#include "flow.h"
// high‐performance push‐relabel
extern "C" {
#include "external/exactDSP-cpp/hi_pr.h"
}

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <numeric>
#include <vector>
using namespace std;

// “infinite” capacity
static const int INF = (int)1e9;
// scaling factor
static const int ACCURACY = 1000;
// max number of density‐improvement rounds
static const int MAX_ROUNDS = 100;
// convergence tolerance
static const double EPS = 1e-16;

// add_arc: add directed arc u->v with given capacity, plus reverse arc v->u of zero capacity
void add_arc(int u, int v, int capacity,
             node *nodes, arc *arcs, cType *cap, long *cur_arc)
{
    arcs[cur_arc[u]].head = nodes + v;
    arcs[cur_arc[u]].rev = arcs + cur_arc[v];
    cap[cur_arc[u]] = capacity;
    ++cur_arc[u];

    arcs[cur_arc[v]].head = nodes + u;
    arcs[cur_arc[v]].rev = arcs + (cur_arc[u] - 1);
    cap[cur_arc[v]] = 0;
    ++cur_arc[v];
}

// nontrivial: returns true if more than one node remains on the source side of the cut
bool nontrivial(int n_nodes, node *nodes)
{
    int res = 0;
    for (int i = 0; i < n_nodes; ++i)
        if (nodes[i].d < n_nodes)
            ++res;
    return res > 1;
}

extern "C" int run_flow(const Graph &G,
                        Stats &stats,
                        const vector<int> &penalty)
{
    int n = G.get_num_nodes();
    int m = G.get_num_edges();
    if (n <= 0 || m <= 0)
    {
        cerr << "Empty graph\n";
        return 1;
    }

    // build undirected edge list (u<v)
    auto adj = G.get_bidirectional_list();
    vector<pair<int, int>> edges;
    edges.reserve(m);
    for (int u = 0; u < n; ++u)
        for (int v : adj[u])
            if (u < v)
                edges.emplace_back(u, v);
    // now edges.size() == m

    // fixed source/sink/node counts
    int source = m + n;
    int sink = source + 1;
    int n_nodes = sink + 1; // = m + n + 2

    // initialize λ
    int sumW = accumulate(penalty.begin(), penalty.end(), 0);
    double lambda = (2.0 * m - sumW) / n;
    cout << lambda << "\n";

    // subgraph indicator
    vector<bool> subg(n, true);

    for (int iter = 0; iter < MAX_ROUNDS; ++iter)
    {
        stats.start_timer();
        // rebuild flow network from scratch each round, using current λ
        // 1) compute prefix‐sum degrees
        vector<long> deg(n_nodes, 0);
        deg[source] = m;
        for (int i = 0; i < m; ++i)
            deg[i] = 3;
        for (int u = 0; u < n; ++u)
            deg[m + u] = 1 + (long)adj[u].size();
        deg[sink] = n;
        for (int i = 1; i < n_nodes; ++i)
            deg[i] += deg[i - 1];
        long n_arcs = deg[n_nodes - 1];

        // allocate
        vector<long> cur_arc(n_nodes);
        vector<bool> subg_next(n, false);
        node *nodes_ = new node[n_nodes + 1];
        arc *arcs = new arc[n_arcs];
        cType *cap = new cType[n_arcs];

        // init pointers
        for (int i = 0; i < n_nodes; ++i)
        {
            cur_arc[i] = (i == 0 ? 0 : deg[i - 1]);
            nodes_[i].first = arcs + cur_arc[i];
        }

        // 2) add arcs
        // -- source → edge nodes
        for (int e = 0; e < m; ++e)
            add_arc(source, e, 2 * ACCURACY, nodes_, arcs, cap, cur_arc.data());
        // -- edge → endpoints
        for (int e = 0; e < m; ++e)
        {
            int u = edges[e].first, v = edges[e].second;
            add_arc(e, m + u, INF, nodes_, arcs, cap, cur_arc.data());
            add_arc(e, m + v, INF, nodes_, arcs, cap, cur_arc.data());
        }
        // -- vertices → sink (using λ)
        for (int u = 0; u < n; ++u)
        {
            int C = int(std::round((lambda + penalty[u]) * ACCURACY));
            add_arc(m + u, sink, C, nodes_, arcs, cap, cur_arc.data());
        }

        // 3) run min‐cut
        min_cut(n_nodes,
                n_arcs / 2,
                nodes_, arcs, cap,
                nodes_ + source, nodes_ + sink,
                0);
        

        // 4) extract new subgraph & compute new λ
        int new_numerator = 0, new_n = 0;
        node *j;
        if (nontrivial(n_nodes, nodes_))
        {
            forAllNodes(j)
            {
                if (nNode(j) < m && nodes_[nNode(j)].d >= n_nodes)
                {
                    int e = nNode(j);
                    new_numerator += 2;
                    subg_next[edges[e].first] = true;
                    subg_next[edges[e].second] = true;
                }
            }
            for (int u = 0; u < n; ++u)
            {
                if (subg_next[u])
                {
                    new_numerator -= penalty[u];
                    ++new_n;
                }
            }
        }
        stats.pause_timer();
        stats.push(iter, lambda, 0.0);
        if (new_n == 0)break;
        double new_c = double(new_numerator) / new_n;
        cout << new_c << "\n";

        delete[] nodes_;
        delete[] arcs;
        delete[] cap;

        // converge?
        if (new_c == lambda)
            break;
        lambda = new_c;
        subg.swap(subg_next);
    }
    std::cout << "Final density: " << lambda << "\n";
    return 0;
}