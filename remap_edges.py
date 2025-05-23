import os
import argparse

"""
Redefine node ids for amazon datasets to be in [0, n), then redefine the mapping.
"""

def remap_edges(input_file, output_file):
    edges = []
    ids = set()

    # Read the file, skip the the num_nodes/edges line
    with open(input_file, 'r') as f:
        first_line = f.readline().strip()
        for line in f:
            parts = line.strip().split()
            if len(parts) != 2:
                continue
            u, v = map(int, parts)
            edges.append((u, v))
            ids.add(u)
            ids.add(v)

    # Map IDs to indices
    id_to_index = {id_: idx for idx, id_ in enumerate(sorted(ids))}

    # Remap edges
    indexed_edges = [(id_to_index[u], id_to_index[v]) for u, v in edges]

    # Sort by first index, then second
    indexed_edges.sort()

    # Write to output file
    with open(output_file, 'w') as f:
        f.write(first_line + '\n')  # Include num_nodes / edges
        for u_idx, v_idx in indexed_edges:
            f.write(f"{u_idx}\t{v_idx}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remap edges")
    parser.add_argument("--path", type=str)
    args = parser.parse_args()

    input_file = os.path.abspath(args.path)
    base, ext = os.path.splitext(os.path.basename(input_file))

    output_dir = os.path.join(os.path.dirname(input_file), 'remap')
    os.makedirs(output_dir, exist_ok=True)  # Create the directory if it doesn't exist
    output_file = os.path.join(output_dir, f"{base}{ext}")

    remap_edges(input_file, output_file)