import os
import random


def build_adjacency_list(file_path):
    adjacency_list = {}

    with open(file_path, "r") as f:
        first_line = f.readline()
        num_nodes, num_edges = map(int, first_line.strip().split())

        # Initialize all node entries
        for i in range(num_nodes):
            adjacency_list[i] = []

        for line in f:
            u, v = map(int, line.strip().split())

            # Store edge only once: store v in u's list, but not u in v's
            adjacency_list[u].append(v)

    return adjacency_list


def find_contra_files(directory):
    contra_files = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith("_contra.txt"):
                full_path = os.path.join(root, file)
                contra_files.append(full_path)
    return contra_files


def read_contra_file(filepath):
    info = []
    with open(filepath, "r") as f:
        for line in f:
            stripped = line.strip()
            if stripped:  # Skip empty lines
                info.append([float(x) for x in stripped.split()])
    return [int(x) for x in info[0]], info[-1][
        0
    ]  # Return all but last line as fixed S, and last line as density


if __name__ == "__main__":
    # ALL _contra.txt files in dsg/dataset/ contain densest sets, with their density on the new line.
    contra_dir = os.path.abspath("./dsg/dataset/remap")
    contra_files = find_contra_files(contra_dir)

    for contra_file_name in contra_files:
        dataset_path = contra_file_name.replace("_contra.txt", ".txt")
        save_path = contra_file_name.replace("_contra.txt", "_perturbed.txt")

        fixed_S, density = read_contra_file(contra_file_name)
        fixed_S = set(fixed_S)

        dataset_path = os.path.abspath(dataset_path)
        adj_list = build_adjacency_list(dataset_path)

        # Compute YES load vector
        load_vector = [0] * len(adj_list)
        # The load for each node in S is density
        for node in fixed_S:
            load_vector[node] = density
        # For every (remaining) edge (u, v). If
        # (u, v) is in E(S) then it's accounted for implicitly.
        # if both (u, v) are outside of S, then we apply 0.5 to both.
        # if one of them is in S, then we apply 1 to the one outside of S.
        for node, neighbors in adj_list.items():
            for neighbor in neighbors:
                if node in fixed_S and neighbor in fixed_S:
                    continue
                if node in fixed_S:
                    load_vector[neighbor] += 1
                elif neighbor in fixed_S:
                    load_vector[node] += 1
                else:
                    load_vector[node] += 0.5
                    load_vector[neighbor] += 0.5

        ground_set = set(range(len(load_vector)))
        S_complement = ground_set - fixed_S

        # Write YES to file
        with open(save_path, "w") as f:
            f.write("Densest set: ")
            f.write(" ".join(map(str, fixed_S)))
            f.write("\nDensity: ")
            f.write(f"{density}\n")
            f.write("YES ")
            f.write(f"{-1} {-1} {-1}: ")
            f.write(" ".join(map(str, load_vector)))

        epsilons = [float(1e-1), float(1), 6, 12]
        for eps in epsilons:
            # Pick arbitrary node in S and one in S complement
            random_in_S = random.choice(list(fixed_S))
            random_out_S = random.choice(list(S_complement))
            n = len(load_vector)

            # Perturb the load vector by eps
            # To create a NO instance, we perturb an element in fixed_S by -eps.
            # We pick arbitrary node in V \ S and perturb it by +eps.
            load_vector[random_in_S] -= eps
            load_vector[random_out_S] += eps

            # Write NO to file
            with open(save_path, "a") as f:
                f.write("\nNO ")
                f.write(f"{eps} {random_in_S} {random_out_S}: ")
                f.write(" ".join(map(str, load_vector)))

            # Recover the load vector
            load_vector[random_in_S] -= eps
            load_vector[random_out_S] += eps

        print("Perturbation done for file:", contra_file_name)
        print("Perturbation saved to file:", save_path)
        print("--------------------------------------------------")
