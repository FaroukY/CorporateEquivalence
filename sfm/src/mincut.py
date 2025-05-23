import os
import re
from collections import deque, defaultdict


def build_adjacency_list_weighted(filename):
    adj = defaultdict(dict)
    with open(filename) as f:
        n, m = map(int, f.readline().strip().split())
        adj = {i: {} for i in range(n)}
        # next(f)  # Skip first line
        for line in f:
            parts = line.strip().split()
            if len(parts) != 3:
                continue
            u, v, w = map(int, parts)
            if u>v:
                continue 
            adj[u][v] = w
            adj[v][u] = w  # undirected
    return adj
def bfs(capacity, source, sink, parent):
    visited = set([source])
    queue = deque([source])
    while queue:
        u = queue.popleft()
        for v in capacity[u]:
            if v not in visited and capacity[u][v] > 0:  # Capacity left on the edge
                queue.append(v)
                visited.add(v)
                parent[v] = u
                if v == sink:
                    return True
    return False
def edmonds_karp(adj, s, t):
    capacity = defaultdict(lambda: defaultdict(int))  # residual capacities
    for u in adj:
        for v, w in adj[u].items():
            capacity[u][v] += w  # Adding weight to handle multiple edges
    parent = {}
    max_flow = 0
    while bfs(capacity, s, t, parent):
        # Find the maximum flow through the path found by BFS
        path_flow = float("Inf")
        v = t
        while v != s:
            u = parent[v]
            path_flow = min(path_flow, capacity[u][v])
            v = u
        # Update residual capacities
        v = t
        while v != s:
            u = parent[v]
            capacity[u][v] -= path_flow
            capacity[v][u] += path_flow
            v = u
        max_flow += path_flow
    return max_flow, capacity
def min_st_cut(adj, s, t):
    # Perform Edmonds-Karp to find the max flow
    max_flow, residual_capacity = edmonds_karp(adj, s, t)
    # Now find the reachable nodes from the source 's' in the residual graph
    visited = set()
    queue = deque([s])
    visited.add(s)
    while queue:
        u = queue.popleft()
        for v in residual_capacity[u]:
            if residual_capacity[u][v] > 0 and v not in visited:
                visited.add(v)
                queue.append(v)
    # Nodes in 'visited' are reachable from s, and form one set of the min-cut
    cut_set_1 = visited
    cut_set_2 = set(residual_capacity.keys()) - visited
    return max_flow, cut_set_1, cut_set_2


if __name__ == "__main__":
    dataset_path = os.path.abspath("/projects/bdzc/eharb/equivalence/sfm/src/dataset/BVZ-sawtooth/remap/BVZ-sawtooth0_s0_t1.max")
    adj_list = build_adjacency_list_weighted(dataset_path)

    # Read s and t from the file name if they match regex
    match = re.search(r's(\d+)_t(\d+)', dataset_path)
    if match:
        s = int(match.group(1))
        t = int(match.group(2))
    else:
        # Default values if not found
        s = 0
        t = len(adj_list) - 1  # Assuming the last node is the sink
    print("Using source:", s, "and sink:", t)
    max_flow, cut_set_1, cut_set_2 = min_st_cut(adj_list, s, t)

    print("Minimum s-t Cut Value (Max Flow):", max_flow)
    # print("Partition 1 (Source side of the cut):", sorted(cut_set_1))