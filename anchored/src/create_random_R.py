import os
import random
import argparse
from collections import defaultdict


def parse_adjacency_list(filename):
    adjacency_list = defaultdict(set)

    with open(filename, "r") as f:
        num_nodes, num_edges = map(int, f.readline().strip().split())

        for line in f:
            if not line.strip():
                continue
            try:
                u, v = map(int, line.strip().split())
                if u != v:  # Skip self-loops
                    adjacency_list[u].add(v)
                    adjacency_list[v].add(u)
            except ValueError:
                continue  # Skip malformed lines

    return adjacency_list


def create_R_set(adj, n_seed, n_total):
    seed_nodes = random.sample(list(adj.keys()), n_seed)
    R_set = set(seed_nodes)
    neighbors = set()
    for node in seed_nodes:
        neighbors.update(adj[node])
    neighbors -= R_set  # avoid reselection

    neighbor_list = list(neighbors)
    random.shuffle(neighbor_list)

    for neighbor in neighbor_list:
        R_set.add(neighbor)
        if len(R_set) >= n_total:
            break

    # If still not enough, fill randomly from remaining graph
    if len(R_set) < n_total:
        remaining = list(set(adj.keys()) - R_set)
        needed = n_total - len(R_set)
        R_set.update(random.sample(remaining, min(needed, len(remaining))))

    return R_set


def main():
    parser = argparse.ArgumentParser(description="Generate R set")
    parser.add_argument("--dataset", type=str, required=True)
    parser.add_argument("--nseed", type=int, default=10)
    parser.add_argument("--ntotal", type=int, default=200)
    args_cmd = parser.parse_args()

    filename = os.path.abspath(args_cmd.dataset)
    adjacency_list = parse_adjacency_list(filename)
    R_set = create_R_set(adjacency_list, args_cmd.nseed, args_cmd.ntotal)

    output_filename = f"{os.path.splitext(filename)[0]}_R_ns{args_cmd.nseed}_nt{args_cmd.ntotal}.txt"
    with open(output_filename, "w") as f:
        f.writelines(f"{node}\n" for node in sorted(R_set))

    print(f"R_set saved to {output_filename}")


if __name__ == "__main__":
    main()
