// High-performance push-relabel DSP solver
extern "C" {
#include "external/exactDSP-cpp/hi_pr.h"
}

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstring> 
#include <iostream>
#include <numeric>
#include <vector>
#include "greedypp.h"
#include "pushrelabel.h" 

using namespace std;
const int INF = (int)1e9;

inline char GET_CHAR()
{
    const int BUFSIZE = 1 << 17;
    static char buf[BUFSIZE], *p = buf, *q = buf;
    if (p == q)
    {
        q = buf + fread(buf, 1, BUFSIZE, stdin);
        p = buf;
        if (p == q)
            return EOF;
    }
    return *p++;
}

inline int getInt() {
    int x = 0;
    char c = GET_CHAR(); 
    
    while ((c < '0' || c > '9') && c != EOF) {
        c = GET_CHAR();
    }

    if (c == EOF) return 0; 

    while (c >= '0' && c <= '9') {
        x = x * 10 + (c - '0');
        c = GET_CHAR(); 
    }
    return x;
}

extern "C" int run_pushrelabel(int argc, char **argv, Stats &stats, Graph &graph)
{
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <ACCURACY> [maxIter]" << endl;
        return 1;
    }
    int ACCURACY = atoi(argv[1]);
    if (ACCURACY <= 0) {
         cerr << "Warning: ACCURACY is non-positive or invalid (" << argv[1] << "), using 100000." << endl;
         ACCURACY = 100000;
    }

    int n = getInt(), m = getInt();
    if (n <= 0 && m > 0) {
        cerr << "Error: n <= 0 but m > 0. (n=" << n << ", m=" << m << ")" << endl;
        return 1;
    }
     if (n < 0 || m < 0) {
        cerr << "Error reading n or m." << endl;
        return 1;
    }

    vector<pair<int, int>> edges(m);
    for (int i = 0; i < m; ++i)
    {
        int u = getInt(), v = getInt();
        if (u < 0 || u >= n || v < 0 || v >= n) {
             cerr << "Error reading edge " << i << " or invalid vertex index (" << u << "," << v << " for n=" << n <<")." << endl;
             return 1;
        }
        edges[i] = {u, v};
    }

    // OPTIMIZATION: Use vector<char> for subg instead of vector<bool>
    vector<char> subg(n, 1); // 1 for true, 0 for false

    // double density = (n > 0) ? (static_cast<double>(m) / n) : 0.0;

    GreedyPeelingPP peeling;
    peeling.run(graph, stats, 1); // One iteration suffices for rough estimate (run counts the time taken and adds it to stats)
    double density = peeling.get_density(); //O(1)

    double prev_density = -1.0; 
    int maxIter = (argc > 2 ? atoi(argv[2]) : 100);
    if (maxIter <= 0) maxIter = 100;


    int max_flow_nodes_alloc = n + m + 2; 
    long max_flow_arcs_alloc = 2L * n + 6L * m + 100; 

    vector<long> deg(max_flow_nodes_alloc);
    vector<long> cur(max_flow_nodes_alloc);
    vector<node> nodes(max_flow_nodes_alloc + 1); 
    vector<arc> arcs(max_flow_arcs_alloc);
    vector<cType> cap(max_flow_arcs_alloc); 

    node *nodes_ptr = nodes.data();
    ::sentinelNode = nodes_ptr + max_flow_nodes_alloc; 

    // new2old and edges1 are still pre-reserved and cleared.
    vector<int> new2old;
    new2old.reserve(n); 
    vector<pair<int, int>> edges1;
    edges1.reserve(m); 

    for (int iter = 0; iter < maxIter; ++iter)
    {
        stats.start_timer();
        prev_density = density;

        // OPTIMIZATION: Revert old2new to per-iteration allocation as in original
        vector<int> old2new(n, -1); // If n=0, this is fine (empty vector)

        new2old.clear();

        int n1 = 0; 
        for (int u = 0; u < n; ++u) {
            if (subg[u]) { // subg[u] is 1 if true
                if (n > 0) old2new[u] = n1; // Check n>0 only if really needed; if subg.size() is n, old2new.size() is also n.
                n1++; // n1 is count, old2new index is n1-1
                new2old.push_back(u);
            }
        }
        // Correct indexing for old2new[u]
        n1 = 0;
        new2old.clear(); // Clear again, as n1 calculation was separate above
        for (int u = 0; u < n; ++u) {
            if (subg[u]) {
                old2new[u] = n1; // n1 starts at 0, so it's the correct new index
                n1++;
                new2old.push_back(u);
            }
        }
        
        if (n1 == 0 && n > 0) { // if n=0, n1 will be 0. If n>0 and n1=0, means empty subgraph.
            density = 0.0;
            if (abs(density - prev_density) < 1e-12 && prev_density != -1.0) {
                // Log and break same as below
                stats.pause_timer(); 
                cerr << "Iter " << iter << ": dens=" << density << " (V=" << 0 << " E=" << 0 << ")\n";
                stats.push(iter, density, 0.0); 
                break;
            }
        }

        edges1.clear();
        for (const auto &e : edges) { 
            int u_new = (n > 0 && e.first < n) ? old2new[e.first] : -1; // defensive access
            int v_new = (n > 0 && e.second < n) ? old2new[e.second] : -1; // defensive access
            if (u_new != -1 && v_new != -1) { 
                edges1.emplace_back(u_new, v_new);
            }
        }

        int m1 = edges1.size();
        int SRC = n1 + m1, SNK = SRC + 1, NND = SNK + 1;

        if (NND > max_flow_nodes_alloc) {
             cerr << "Error: NND " << NND << " exceeds allocation " << max_flow_nodes_alloc << endl; return 1;
        }
        // Ensure deg is filled only up to NND
        fill(deg.begin(), deg.begin() + NND, 0L);
        if (NND > 0) { // Only assign if NND is valid
            if (SRC < NND) deg[SRC] = n1;
            if (SNK < NND) deg[SNK] = m1; 
        }

        for (int i = 0; i < n1; ++i) {
            if (i < NND) deg[i] += 1; // Original was deg[i]=1, if NND is small, +=1 could be an issue if not init to 0.
                                     // With fill above, it's fine.
        }
        for (const auto &e1 : edges1) { 
            if (e1.first < NND) deg[e1.first]++;
            if (e1.second < NND) deg[e1.second]++;
        }
        for (int j = 0; j < m1; ++j) {
            if ((n1 + j) < NND) deg[n1 + j] += 3; // Original was deg[n1+j]=3
        }
        // To be absolutely safe and match original logic for deg array assignments
        // (assuming deg was filled with 0s up to NND before this):
        // Re-do the deg assignment exactly as original:
        if (NND > 0) {
            if (SRC < NND) deg[SRC] = n1; else if (SRC==0 && n1==0) deg[SRC]=0; // handle NND=1, SRC=0 case etc.
            if (SNK < NND) deg[SNK] = m1; else if (SNK==0 && m1==0) deg[SNK]=0;
        }
        for (int i = 0; i < n1; ++i) { // These are vertex nodes 0..n1-1
             deg[i] = 1; 
        }
        for (const auto &e1 : edges1) { 
            deg[e1.first]++; // e1.first is a new vertex index < n1
            deg[e1.second]++;// e1.second is a new vertex index < n1
        }
        for (int j = 0; j < m1; ++j) { // These are edge nodes n1..n1+m1-1
            deg[n1 + j] = 3;
        }


        for (int i = 1; i < NND; ++i)
            deg[i] += deg[i - 1];
        long tot_structural_arcs = (NND > 0) ? deg[NND - 1] : 0;

        if (tot_structural_arcs > max_flow_arcs_alloc) {
            cerr << "Error: tot_structural_arcs " << tot_structural_arcs << " exceeds allocation " << max_flow_arcs_alloc << endl; return 1;
        }

        if (NND > 0) cur[0] = 0;
        for (int i = 1; i < NND; ++i)
            cur[i] = deg[i - 1];
        
        for (int i = 0; i < NND; ++i)
            nodes_ptr[i].first = arcs.data() + cur[i];

        auto add_arc = [&](int u_arc, int v_arc, cType capacity_val) { 
            long pu = cur[u_arc]++; 
            long pv = cur[v_arc]++; 
            
            arcs[pu].head = &nodes_ptr[v_arc];
            arcs[pu].rev = &arcs[pv];
            cap[pu] = capacity_val;

            arcs[pv].head = &nodes_ptr[u_arc];
            arcs[pv].rev = &arcs[pu];
            cap[pv] = 0; 
        };
        
        cType src_u_cap_val = static_cast<cType>(density * ACCURACY);
        if (src_u_cap_val < 0) src_u_cap_val = 0;


        for (int u_new = 0; u_new < n1; ++u_new)
            add_arc(SRC, u_new, src_u_cap_val);
        
        for (int j = 0; j < m1; ++j) {
            const auto &e1 = edges1[j]; 
            add_arc(e1.first, n1 + j, INF);  
            add_arc(e1.second, n1 + j, INF); 
        }
        for (int j = 0; j < m1; ++j)
            add_arc(n1 + j, SNK, static_cast<cType>(ACCURACY));

        if (NND > 0 && tot_structural_arcs > 0 && SRC < NND && SNK < NND && SRC != SNK) { 
             min_cut(NND, tot_structural_arcs / 2, nodes_ptr, arcs.data(), cap.data(), 
                    &nodes_ptr[SRC], &nodes_ptr[SNK], 0);
        }
        
        if (n > 0) fill(subg.begin(), subg.end(), (char)0); // Fill with char 0
        long src_arc_base_offset = (NND > 0 && SRC < NND && nodes_ptr[SRC].first != nullptr) ? (nodes_ptr[SRC].first - arcs.data()) : -1;
        
        int vcount = 0;
        if (src_arc_base_offset != -1) { 
            for (int idx = 0; idx < n1; ++idx) { 
                // Check node_ptr[idx].d is initialized by min_cut if NND is small.
                // min_cut should initialize d for all NND nodes.
                if (nodes_ptr[idx].d < NND && cap[src_arc_base_offset + idx] > 0) {
                    if (n > 0) subg[new2old[idx]] = 1; // 1 for true
                    vcount++;
                }
            }
        }
        
        int ecount = 0; 
        for(const auto& edge_pair_new_indices : edges1) { 
            if (n > 0 && subg[new2old[edge_pair_new_indices.first]] && subg[new2old[edge_pair_new_indices.second]]) {
                ecount++;
            }
        }

        density = (vcount > 0) ? (static_cast<double>(ecount) / vcount) : 0.0;

        stats.pause_timer(); 
        cerr << "Iter " << iter << ": dens=" << density << " (V=" << vcount << " E=" << ecount << ")\n";
        stats.push(iter, density, 0.0); 

        if (prev_density >= 0.0 && abs(density - prev_density) < 1e-12) 
            break;
        // Special case: if n=0, density is 0. if prev_density was also 0 (or -1 initially), break.
        if (n == 0 && density == 0.0 && (prev_density == 0.0 || prev_density == -1.0)) break;
    }

    return 0;
}