import os
from collections import defaultdict
import random
import argparse

def build_adjacency_list_weighted(filename):
    adj = defaultdict(dict)
    edges = []
    with open(filename) as f:
        next(f)  # Skip first line
        for line in f:
            parts = line.strip().split()
            if len(parts) != 3:
                continue
            u, v, w = map(int, parts)
            if u>v:
                continue
            adj[u][v] = w
            adj[v][u] = w  # undirected
            edges.append((u, v, w))
    return adj, edges


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate R set")
    parser.add_argument("--dataset", type=str, default="./sfm/src/dataset/m50.txt")
    parser.add_argument("--n_nodes_to_add", type=int, default=100) # vary from 100 to 100_000?
    parser.add_argument("--n_edges_to_add", type=int, default=-1)
    args_cmd = parser.parse_args()

    dataset_path = os.path.abspath(args_cmd.dataset)
    adj_list, edges = build_adjacency_list_weighted(dataset_path)
    s = 0
    t = len(adj_list) - 1  # Assuming the last node is the sink

    max_edge_weight = max(max(w.values()) for w in adj_list.values())

    n_nodes = args_cmd.n_nodes_to_add + len(adj_list)

    edges_added = 0
    edges_to_add = args_cmd.n_edges_to_add if args_cmd.n_edges_to_add != -1 else n_nodes * 2

    while edges_added < edges_to_add:
        # Randomly select node from adj_list (in graph)
        u = random.choice(list(adj_list.keys()))
        # Randomly select node from new nodes (not in graph)
        v = random.randint(0, n_nodes - 1)
        # Weight for the edge
        w = random.randint(1, max_edge_weight)

        if u == v or adj_list[u].get(v) is not None or adj_list[v].get(u) is not None:
            # Edge already exists, skip this iteration
            continue
        adj_list[u][v] = w
        adj_list[v][u] = w  # undirected
        edges.append((u, v, w))
        edges_added += 1

    # Remap every node to a new index, starting from 0
    # and remap every edge to these new indices
    node_mapping = {}
    new_index = 0
    nodes = list(adj_list.keys())
    for node in nodes:
        node_mapping[node] = new_index
        new_index += 1
    new_edges = []
    for u, v, w in edges:
        new_u = node_mapping[u]
        new_v = node_mapping[v]
        new_edges.append((new_u, new_v, w))

    # Save the modified adjacency list to a new file with s, t as the first line
    output_file = os.path.abspath(f"./sfm/src/dataset/m50_random_s{s}_t{t}.txt")
    num_nodes = len(adj_list)

    new_edges.sort(key=lambda x: (x[0], x[1]))

    with open(output_file, "w") as f:
        f.write(f"{num_nodes} {len(new_edges)}\n")
        for u, v, w in new_edges:
            f.write(f"{u} {v} {w}\n")
    print(f"Random graph saved to {output_file}")
    print(f"Number of nodes: {num_nodes}")
    print(f"Number of edges: {len(new_edges)}")

