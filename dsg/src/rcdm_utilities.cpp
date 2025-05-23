#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <numeric>
#include <queue>
#include <utility>
#include <vector>
using namespace std;
#include "rcdm_utilities.h"
/*
This uses fractional peeling...
*/
double current_density(vector<vector<pair<int, int>>> &Adj, vector<double> &b, vector<double> &x, long long n, long long m, long long t)
{
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;

    for (int i = 0; i < n; ++i)
        pq.push({b[i] + Adj[i].size(), i});
    vector<bool> deleted(n, false);
    vector<double> degree(n, 0);
    for (int i = 0; i < n; ++i)
        degree[i] = b[i] + Adj[i].size();
    long long N = n;
    long long M = m;
    double density = (double)(1.0 * M / N);
    while (!pq.empty())
    {
        auto [d, i] = pq.top();
        pq.pop();
        if (deleted[i])
            continue;
        for (auto &[j, e_idx] : Adj[i])
        {
            int sister_idx = e_idx % 2 ? e_idx - 1 : e_idx + 1;
            if (!deleted[j])
            {
                pq.push({degree[j] - t * x[sister_idx] - 1, j});
                degree[j] -= (t * x[sister_idx] + 1);
                M--;
            }
        }
        N--;
        if (N > 0)
        {
            density = max(density, 1.0 * M / N);
        }
        deleted[i] = true;
    }
    return density;
}

double current_density_sorting(vector<vector<pair<int, int>>> &Adj, vector<double> &b, long long n, long long m)
{
    vector<int> indices(n);
    vector<bool> deleted(n, false);
    iota(indices.begin(), indices.end(), 0);

    sort(indices.begin(), indices.end(), [&b](const int &i, const int &j) -> bool { return b[i] < b[j]; });
    double density = 1.0 * m / n;
    long long N = n, M = m;
    for (auto &i : indices)
    {
        for (auto &[j, e_idx] : Adj[i])
        {
            if (!deleted[j])
                M--;
        }
        --N;
        if (N)
            density = max(density, 1.0 * M / N);
        deleted[i] = true;
    }
    return density;
}
vector<double> get_initial(vector<vector<pair<int, int>>> &Adj, long long n, long long m)
{
    vector<double> x(2 * m, 0.0);

    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    for (int i = 0; i < n; ++i)
        pq.push({Adj[i].size(), i});
    vector<bool> deleted(n, false);
    vector<int> degree(n, 0);
    for (int i = 0; i < n; ++i)
        degree[i] = Adj[i].size();
    long long N = n;
    long long M = m;
    double density = 1.0 * M / N;
    while (!pq.empty())
    {
        auto [d, i] = pq.top();
        pq.pop();
        if (deleted[i])
            continue;

        for (auto &[j, e_idx] : Adj[i])
        {
            if (!deleted[j])
            {
                if (i != j)
                {
                    pq.push({--degree[j], j});
                    M--;
                }
                x[e_idx] = 1;
            }
        }
        N--;
        if (N > 0)
        {
            density = max(density, 1.0 * M / N);
        }
        deleted[i] = true;
    }
    // assert (abs(accumulate(x.begin(), x.end(), 0)-m)<=1e-5 );
    return x;
}