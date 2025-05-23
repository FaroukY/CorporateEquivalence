def convert_metis_to_edge_list(input_file, output_file):
    edges = set()  # Use a set to avoid duplicates (for undirected edges)
    used_vertices = set()

    with open(input_file, 'r') as f:
        lines = f.readlines()

    # Remove comment lines (starting with %)
    lines = [line.strip() for line in lines if line.strip() and not line.strip().startswith('%')]

    # Read first line: n and m (we'll recompute n later after removing isolated vertices)
    header = lines[0].split()
    original_n = int(header[0])

    # Read adjacency lists
    for u, line in enumerate(lines[1:], start=0):  # start=0 for 0-based index
        neighbors = list(map(int, line.split()))
        for v in neighbors:
            v0 = v - 1  # Convert to 0-based index
            if u <= v0:  # Avoid duplicating undirected edges
                edges.add((u, v0))
                used_vertices.add(u)
                used_vertices.add(v0)

    # Build remapping from old vertex ID to new ID (0 to n'-1)
    used_vertices = sorted(used_vertices)
    vertex_id_map = {old_id: new_id for new_id, old_id in enumerate(used_vertices)}

    # Re-map edges
    remapped_edges = {(vertex_id_map[u], vertex_id_map[v]) for u, v in edges}

    # Write to output file
    n_prime = len(used_vertices)
    with open(output_file, 'w') as f:
        f.write(f"{n_prime} {len(remapped_edges)}\n")
        for u, v in sorted(remapped_edges):
            f.write(f"{u} {v} 1\n")

# Example usage:
convert_metis_to_edge_list('kron_g500-simple-logn19.graph', 'kron_g500_logn19_cleaned.txt')
