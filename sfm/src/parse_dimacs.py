import argparse
import os

def parse_dimacs_to_undirected(input_path, output_path):
    nodes = set()
    edge_map = {}
    s, t = None, None

    with open(input_path, "r") as infile:
        for line in infile:
            line = line.strip()
            if line.startswith("a"):
                # arc line: a u v capacity [optional flow]
                parts = line.split()
                if len(parts) < 4:
                    raise ValueError(f"Unexpected arc line format (too short): {line}")
                _, u_str, v_str, capacity_str = parts[:4]  # Only take first 4 fields
                u = int(u_str) - 1  # Convert to 0-based index
                v = int(v_str) - 1
                capacity = int(capacity_str)

                # Normalize undirected edge as (min, max)
                edge = (min(u, v), max(u, v))

                if edge in edge_map:
                    edge_map[edge] += capacity
                else:
                    edge_map[edge] = capacity

                nodes.add(u)
                nodes.add(v)

            elif line.startswith("n"):
                # source/sink line: n <node_id> s|t
                parts = line.split()
                if len(parts) != 3:
                    raise ValueError(f"Unexpected node line format: {line}")
                _, node_str, node_type = parts
                node = int(node_str) - 1  # Convert to 0-based index

                if node_type == 's':
                    s = node
                elif node_type == 't':
                    t = node

    if s is None or t is None:
        raise ValueError("Source (s) or sink (t) node not found in input.")

    num_nodes = max(nodes) + 1
    edges = sorted(edge_map.items(), key=lambda x: (x[0][0], x[0][1]))

    with open(output_path, "w") as outfile:
        outfile.write(f"{num_nodes} {len(edges)}\n")
        for (u, v), cap in edges:
            outfile.write(f"{u} {v} {cap}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert DIMACS arcs to undirected edges.")
    parser.add_argument("--path", type=str, required=True, help="Path to input DIMACS .max file")
    args = parser.parse_args()

    input_file = os.path.abspath(args.path)
    base, ext = os.path.splitext(os.path.basename(input_file))

    # First pass to extract s, t node ids
    s, t = None, None
    with open(input_file, "r") as infile:
        for line in infile:
            if line.startswith("n"):
                _, node_str, node_type = line.strip().split()
                node = int(node_str) - 1  # Convert to 0-based index
                if node_type == 's':
                    s = node
                elif node_type == 't':
                    t = node
    if s is None or t is None:
        raise ValueError("Source (s) or sink (t) node not found in input.")

    # Build output filename with _s{source}_t{sink}
    output_dir = os.path.join(os.path.dirname(input_file), "remap")
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{base}_s{s}_t{t}{ext}")

    # Run main conversion
    parse_dimacs_to_undirected(input_file, output_file)

    print(f"Converted file saved to {output_file}")
