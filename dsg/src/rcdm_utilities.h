#pragma once
#include <vector>

double current_density(std::vector<std::vector<std::pair<int, int>>> &Adj, std::vector<double> &b, std::vector<double> &x, long long n, long long m, long long t = 1);

double current_density_sorting(std::vector<std::vector<std::pair<int, int>>> &Adj, std::vector<double> &b, long long n, long long m);

std::vector<double> get_initial(std::vector<std::vector<std::pair<int, int>>> &Adj, long long n, long long m);