#include "greedypp.h"

#include "logger.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

using namespace std;

namespace {
    inline char GET_CHAR()
    {
        const int maxn = 131072;
        static char buf[maxn], *p1 = buf, *p2 = buf;
        return p1 == p2 && (p2 = (p1 = buf) + fread(buf, 1, maxn, stdin), p1 == p2) ? EOF : *p1++;
    }

    inline int getInt()
    {
        int res(0);
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

    struct Edge
    {
        int y, next;
    };
    struct Node
    {
        int deg, next, prev, idx;
        void clear()
        {
            deg = next = prev = 0;
            idx = -1;
        }
    };

    inline void linklists(Node *lists, int x, int y)
    {
        if (y == 0)
            return;
        lists[x].next = y;
        lists[y].prev = x;
    }

    inline void linknodes(int *nxt, int *prv, int x, int y)
    {
        if (y == -1)
            return;
        nxt[x] = y;
        prv[y] = x;
    }

    inline void eraselist(Node *lists, int x)
    {
        lists[lists[x].prev].next = lists[x].next;
        if (lists[x].next != 0)
            lists[lists[x].next].prev = lists[x].prev;
    }

    inline void erasenode(Node *lists, int *nxt, int *prv, int *itr, int x)
    {
        if (prv[x] == -1)
        {
            lists[itr[x]].idx = nxt[x];
        }
        if (prv[x] != -1)
            nxt[prv[x]] = nxt[x];
        if (nxt[x] != -1)
            prv[nxt[x]] = prv[x];
    }
} // namespace

void GreedyPeelingPP::run(const Graph &graph, Stats &stats, int iterations)
{
    // Original code with minimal structure change
    int n = graph.get_num_nodes();
    int m = graph.get_num_edges();

    vector<vector<pair<int, int>>> Adj = graph.get_adjacency_list();

    vector<Edge> edges(m * 2 + 10);
    vector<int> idx(n, 0);
    vector<int> init_deg(n, 0);

    int l = 0;
    auto build = [&](int x, int y) {
        edges[++l].next = idx[x];
        edges[l].y = y;
        idx[x] = l;
    };

    for (int u = 0; u < n; ++u)
    {
        for (auto &[v, _] : Adj[u])
        {
            build(u, v);
            init_deg[u]++;
        }
    }

    vector<Node> lists(n + 2 * m + 10);

    w = vector<int>(n, 0);
    vector<int> deg(n), pos(n), itr(n), prv(n), nxt(n);
    vector<pair<int, int>> deg_sorted(n);

    mm_density = 0;
    vector<int> m_ans;

    for (int tt = 0; tt < iterations; ++tt)
    {
        stats.start_timer();

        for (int i = 0; i < n; ++i)
        {
            nxt[i] = prv[i] = -1;
            pos[i] = 0;
            deg[i] = w[i] + init_deg[i];
            deg_sorted[i] = {deg[i], i};
        }

        sort(deg_sorted.begin(), deg_sorted.end());
        int n_list = 0;

        for (int i = 0; i < n; ++i)
        {
            int v = deg_sorted[i].second;
            if (n_list == 0 || lists[n_list].deg != deg_sorted[i].first)
            {
                ++n_list;
                lists[n_list].clear();
                linklists(lists.data(), n_list - 1, n_list);
                lists[n_list].deg = deg_sorted[i].first;
            }
            linknodes(nxt.data(), prv.data(), v, lists[n_list].idx);
            lists[n_list].idx = v;
            itr[v] = n_list;
        }

        double max_density = (double)m / n;
        int cur_m = m, cur_n = n;
        vector<int> ans;
        int max_size = 0;

        while (lists[0].next)
        {
            int i = lists[0].next;
            int k = lists[i].idx;

            if (nxt[k] == -1)
                eraselist(lists.data(), i);
            else
                erasenode(lists.data(), nxt.data(), prv.data(), itr.data(), k);

            pos[k] = -1;
            w[k] = deg[k];
            cur_n -= 1;
            ans.push_back(k);

            for (int p = idx[k]; p; p = edges[p].next)
            {
                int j = edges[p].y;
                if (pos[j] == -1)
                    continue;
                cur_m -= 1;

                int i1 = itr[j];
                erasenode(lists.data(), nxt.data(), prv.data(), itr.data(), j);
                int i_prev = lists[i1].prev;

                if (lists[i1].idx == -1)
                    eraselist(lists.data(), i1);
                deg[j]--;

                prv[j] = nxt[j] = -1;
                if (i_prev == 0 || lists[i_prev].deg != deg[j])
                {
                    ++n_list;
                    lists[n_list].clear();
                    itr[j] = n_list;
                    int next_in_list = lists[i_prev].next;
                    lists[n_list].deg = deg[j];
                    lists[n_list].idx = j;
                    linklists(lists.data(), i_prev, n_list);
                    if (next_in_list)
                        linklists(lists.data(), n_list, next_in_list);
                } else
                {
                    linknodes(nxt.data(), prv.data(), j, lists[i_prev].idx);
                    lists[i_prev].idx = j;
                    itr[j] = i_prev;
                }
            }

            if (cur_n == 0)
                continue;
            if (max_density < (double)cur_m / cur_n)
                max_size = ans.size();

            max_density = max(max_density, (double)cur_m / cur_n);
        }

        reverse(ans.begin(), ans.end());
        ans.resize(n - max_size);

        if (max_density > mm_density)
        {
            m_ans = ans;
            mm_density = max_density;
        }

        stats.pause_timer();

        double loadnorm = 0.0;
        for (int i = 0; i < n; ++i)
            loadnorm += w[i] * w[i] / ((tt + 1) * (tt + 1));
        log("Iteration ", tt + 1, ": ", mm_density);
        stats.push(tt, mm_density, loadnorm);
    }
}