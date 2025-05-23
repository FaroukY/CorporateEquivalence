# Equivalence Project

## About 

This repository contains the codebase for the paper "Corporate Needs You to Find the Difference: Revisiting Submodular and Supermodular Ratio Optimization Problems." All experiments presented in the main paper can be reproduced using the code and instructions provided in this README.

## Accessing Full and Cleaned Datasets. 

This repository includes only one representative dataset per problem to help users quickly run the codebase. To access the full set of datasets used in the paper, please download them from the following (anonymous) link:

https://drive.google.com/file/d/1kuAm4r_0PfRasAfdwudMldvxd368n7yC/view?usp=sharing. 

After downloading, extract the contents and move the datasets into their respective directories in this code-base (e.g., `dsg/dataset/remap/`). This will enable you to reproduce all experiments presented in the paper.

## Required Libraries. 

To run the experiments in this repository, ensure the following libraries are installed:

1. Eigen: The C++ linear algebra library (version 3.3 or higher is required).
2. Boost: The C++ boost library. We used boost_1_81_0. 
2. TBB (Threading Building Blocks): Required for efficient parallel computation.
3. Gurobi (for HNSN experiments in the usss directory only): Make sure to set the Gurobi license file environment variable: `export GRB_LICENSE_FILE=/sw/external/gurobi/gurobi.lic`

You have the choice of downloading the above on your system, or building from source for Eigen and Boost. To build Eigen from source, include the source under the `/external/eigen` path. By default, we fetch and link Boost from source; to disable this behaviour disable `USE_FETCHCONTENT_FOR_BOOST` in the make CMake file by changing `ON` to `OFF`.

If you encounter an error such as `libtbb.so.12: cannot open shared object file: No such file or directory` when running the experiment driver, set the LD_LIBRARY_PATH environment variable to include the path to your TBB installation (adjust the path as needed based on where you downloaded TBB):

```bash
export LD_LIBRARY_PATH=/sw/spack/deltas11-2023-03/apps/linux-rhel8-zen/gcc-8.5.0/intel-tbb-2021.9.0-m5eh674/lib64:$LD_LIBRARY_PATH
```

These dependencies must be properly set up to compile and run the full set of experiments.

## Build Instructions

To compile the entire project:

```bash
 cmake --preset=dev .
 cmake --build --preset=dev  
```

To compile a specific target (say dsg), specify the `<project>_driver` as well:
```bash
 cmake --preset=dev .
 cmake --build --preset=dev --target dsg_driver    
```

## Slurm Generation

The experiments for a project are defined in an accordingly-named yaml file within `configs/`. Each
config specifies a set of default resource parameters, alongside algorithms and datasets whose cartesian product defines the set of experiments. To generate the corresponding slurms, run the generation script pointing to the desired config. For example, `dsg` slurms are generated using:

```bash
python3 ./configs/generate_slurms.py --config=./configs/dsg.yaml --slurm_dir=slurms/dsg
```

You can then submit them to your cluster using:

```bash
for file in slurms/dsg/*.slurm; do sbatch "$file"; done
```

To modify the template, you can modify the slurm template in `./configs/generate_slurms.py`. 

## Suppressing Logs

Compile with `ENABLE_LOGGING=0` in the relevant `CMakeLists` to disable logs; use `1` to enable them.

---

## Densest Subgraph (DSG)

### Overview

Implemented algorithms for finding densest subgraphs:

1. **FISTA** – Convex optimization, in `dsg/src/fista.h`
2. **Frank-Wolfe** – Convex optimization, in `dsg/src/frankwolfe.h`
3. **Minimum Norm Point (MNP)** – Convex optimization, in `dsg/src/mnp.h`
4. **RCDM** -- Convex optimization, in `dsg/src/rcdm_permutation.h`
5. **Greedy++** -- Combinatorial, in `dsg/src/greedypp.h`
6. **Incremental Flow** – Flow-based, in `dsg/src/incremental.h`
7. **Push-Relabel (NEW ALGORITHM)** – Flow-based, in `dsg/src/pushrelabel.h`

### Datasets

Available in the `dsg/dataset/remap` directory:
- `close-cliques.txt`
- `roadNet-PA.txt`
- `cit-Patents.txt`
- `com-amazon.ungraph.txt`
- `com-dblp.ungraph.txt`
- `com-orkut.ungraph.txt`
- `roadNet-CA.txt`

Expected input format for graphs:
```
n m
u1 v1
u2 v2
...
um vm
```
Where `n` = number of nodes, `m` = number of edges, and `0 ≤ ui, vi < n`. No self-loops.

### Running the Experiments

All algorithms are integrated via `dsg/src/main.cpp`. While in the `root` directory, run:

```bash
./build/dev/dsg/dsg_driver <PATH_TO_DATASET.txt> <ALGORITHM_NAME> <OUTPUT_PATH.txt> <OPTIONAL_ARGS>
```

#### Examples

```bash
./build/dev/dsg/dsg_driver ./dsg/dataset/remap/close-cliques.txt fista ./output/dsg/out.txt 5 #5 iterations
./build/dev/dsg/dsg_driver ./dsg/dataset/remap/close-cliques.txt frankwolfe ./output/dsg/out.txt 5
./build/dev/dsg/dsg_driver ./dsg/dataset/remap/close-cliques.txt mnp ./output/dsg/out.txt 5
./build/dev/dsg/dsg_driver ./dsg/dataset/remap/close-cliques.txt rcdm ./output/dsg/out.txt 5
./build/dev/dsg/dsg_driver ./dsg/dataset/remap/close-cliques.txt greedypp ./output/dsg/out.txt 5
./build/dev/dsg/dsg_driver ./dsg/dataset/remap/close-cliques.txt incremental ./output/dsg/out.txt
./build/dev/dsg/dsg_driver ./dsg/dataset/remap/close-cliques.txt pushrelabel ./output/dsg/out.txt
```

For the experiments above, the file ./output/dsg/out.txt contains one row per iteration of the algorithm. Each entry includes, in order: the iteration number, the density for that iteration, the norm of the load vector (for convex optimization-based methods; otherwise, this will be 0), and the time taken for that iteration (in nanoseconds, not cumulative).

---

### Anchored Densest Subgraph (ADSG)

### Overview

Implemented algorithms for finding anchored densest subgraphs:

2. **Frank-Wolfe** – Convex optimization, in `anchored/src/frankwolfe.h`
3. **Minimum Norm Point (MNP)** – Convex optimization, in `anchored/src/mnp.h`
5. **SuperGreedy++** -- Combinatorial, in `anchored/src/supergreedy.h`
7. **Flow Algorithm** – Flow-based, in `anchored/src/flow.h`

### Datasets

We reuse the datasets from the $p$-means densest subgraph problem (DSS). See details below. Specifically, the datasets are in `dsg/dataset/remap/`. 

To create an R_set for a given DSG dataset, point to it and run

```bash
python3 ./anchored/src/create_random_R.py --dataset ./dsg/dataset/remap/close-cliques.txt --nseed 10 --ntotal 201
``` 

which uses a seed of 10 and grows the set to |R| = 201.

The R set is saved to DSG datasets and suffixed by the seed and total. An ADSG experiment can then be run with the original DSG dataset by pointing to the corresponding anchor set to use:

### Running the Experiments

All algorithms are integrated via `anchored/src/main.cpp`. While in the `root` directory:

```bash
./build/dev/anchored/anchored_driver <PATH_TO_DATASET.txt> <ALGORITHM_NAME> <OUTPUT_PATH.txt> <OPTIONAL_ARGS>
```

### Examples
```bash
./build/dev/anchored/anchored_driver ./dsg/dataset/remap/close-cliques_R_ns10_nt201.txt frankwolfe ./output/adsg/out.txt 100
./build/dev/anchored/anchored_driver ./dsg/dataset/remap/close-cliques_R_ns10_nt201.txt supergreedy ./output/adsg/out.txt 100
./build/dev/anchored/anchored_driver ./dsg/dataset/remap/close-cliques_R_ns10_nt201.txt mnp ./output/adsg/out.txt 100
./build/dev/anchored/anchored_driver ./dsg/dataset/remap/close-cliques_R_ns10_nt201.txt flow ./output/adsg/out.txt
```

For the experiments above, the file `./output/adsg/out.txt` contains one row per iteration of each  algorithm. Each entry includes, in order: the iteration number, the density for that iteration, the norm of the load vector (for convex optimization-based methods; otherwise for flow, this will be 0), and the time taken for that iteration (in nanoseconds, not cumulative).

---

## Densest Supermodular Set (DSS)

### Overview 

Objective: Maximize the generalized density `f(S)/|S|`, where `f(S)` is a normalized supermodular function (`f(∅) = 0`).

A primary instance is:

$$
f(S) = \sum_{u ∈ S} \deg_S(u)^p,~ p > 1
$$

This generalizes densest subgraph computation. Flow-based solutions don't apply here.

We implement the following algorithms:

- Frank-Wolfe, Implemented in `dss/src/frankwolfe.h`
- SuperGreedy++, Implemented in `dss/src/supergreedy.h`
- Minimum Norm Point (MNP), Implemented in `dss/src/mnp.h`

### Datasets

We use a subset of the datasets we used for DSG, available in the `dsg/dataset/remap` directory:
- `close-cliques.txt`
- `roadNet-PA.txt`
- `com-amazon.ungraph.txt`
- `com-dblp.ungraph.txt`
- `roadNet-CA.txt`

We use $p=1.1, 1.25, 1.5, 1.75$. 

### Running the Experiments

All algorithms are integrated via `dss/src/main.cpp`. While in the `root` directory:

```bash
./build/dev/dss/dss_driver <PATH_TO_DATASET.txt> <ALGORITHM_NAME> <OUTPUT_PATH.txt> <NUMBER_OF_ITERATIONS> <P_VALUE>
```

### Examples
```bash
./build/dev/dss/dss_driver ./dsg/dataset/remap/close-cliques.txt frankwolfe ./output/dss/out.txt 100 1.1
./build/dev/dss/dss_driver ./dsg/dataset/remap/close-cliques.txt supergreedy ./output/dss/out.txt 100 1.25
./build/dev/dss/dss_driver ./dsg/dataset/remap/close-cliques.txt mnp ./output/adsg/dss.txt 100 1.75
```

For the experiments above, the file `./output/dss/out.txt` contains one row per iteration of each  algorithm. Each entry includes, in order: the iteration number, the density for that iteration, the norm of the load vector, and the time taken for that iteration (in nanoseconds, not cumulative).

---

## Unrestricted Sparsest Submodular Set (USSS) (HNSN Problem)

### Overview

For USSS, we select the problem of Heavy Nodes in a Small Neighborhood (HNSN) from the paper
Heavy Nodes in a Small Neighborhood: Exact and Peeling Algorithms with Applications 

The problem is given a bipartite graph as input (L, R), with weights on right vertices
represented as $w(r), r\in R$. Goal is to select a set $S \subseteq R$ that minimizes 
$|N(S)|/w(S)$. $f(S)=|N(S)|$ is obviously submodular, and $w(S)$ is modular. 

Implemented algorithms for HNSN (with the letters inside paranthesis indicating algorithm name):

1. **SuperGreedy++ Iterative Peeling (IP)** – Convex optimization, in `usss/src/supergreedy.h`
2. **Frank-Wolfe (FW)** – Convex optimization, in `usss/src/frankwolfe.h`
3. **Minimum Norm Point (MNP)** – Convex optimization, in `usss/src/mnp.h`
4. **Linear Programming (LP)** - Linear Programming, in `usss/src/main.cpp` Implemented by Ling et al. 
5. **ContractDecompose (CD)** - Combinatorial, in `usss/src/main.cpp` Implemented by Ling et al., due to Ling et al. 
6. **Greedy Approximation (GAR)** - Combinatorial, in `usss/src/main.cpp` Implemented by Ling et al. 
7. **Greedy (GR)** - Combinatorial, in `usss/src/main.cpp` Implemented by Ling et al. 
8. **Fast Greedy (FGR)** - Combinatorial, in `usss/src/main.cpp` Implemented by Ling et al., due to Ling et al. 
9. **GreedRatio (GRR)** - Combinatorial, in `usss/src/main.cpp` Implemented by Ling et al. due to W. Bai et al. 
10. **NEW Flow Based Algorithm (FLOW)** - Combinatorial, in `usss/src/flow.h`. 


### Datasets

Available in the `usss/datasets/` directory:
- `graph_NBA_updated`
- `graph_foodmart_updated`
- `graph_ecommerce_updated`
- `graph_liquor_updated`
- `graph_yoochoose_updated`
- `graph_fruithut_updated`
- `graph_connectious_updated`
- `graph_notredame_updated`
- `graph_ACM_updated`
- `graph_imdb_updated`
- `graph_kosarak_updated`
- `graph_digg_updated`

Expected input format for bipartite graphs:

**IMPORTANT: Right vertices must be labeled 0 to numRight-1, and Left vertices numRight to numRight+numLeft**

```
numLeft
numRight
weights_{1}
...
weights_{numRight}
u1 v1
u2 v2
...
um vm
```
Where `numLeft` = number of nodes on left, `numRight` = number of nodes on right, `m` = number of edges. 

### Running the Experiments

All algorithms are integrated via `usss/src/main.cpp`. While in the `root` directory:

```bash
./build/dev/usss/usss_driver <DATASET_PATH> <NAME_OF_ALGORITHM> <OUTPUT_PATH.txt> <OPTIONAL_ITERATIONS>
```

#### Examples

```bash
./build/dev/usss/usss_driver usss/dataset/graph_ecommerce_updated IP output/hnsn/out.txt 100
./build/dev/usss/usss_driver usss/dataset/graph_ecommerce_updated FW output/hnsn/out.txt 100
./build/dev/usss/usss_driver usss/dataset/graph_ecommerce_updated MNP output/hnsn/out.txt 100
export GRB_LICENSE_FILE=/sw/external/gurobi/gurobi.lic && ./build/dev/usss/usss_driver usss/dataset/graph_ecommerce_updated LP output/hnsn/out.txt
export GRB_LICENSE_FILE=/sw/external/gurobi/gurobi.lic && ./build/dev/usss/usss_driver usss/dataset/graph_ecommerce_updated CD output/hnsn/out.txt
./build/dev/usss/usss_driver usss/dataset/graph_ecommerce_updated GAR output/hnsn/out.txt
./build/dev/usss/usss_driver usss/dataset/graph_ecommerce_updated GR output/hnsn/out.txt
./build/dev/usss/usss_driver usss/dataset/graph_ecommerce_updated FGR output/hnsn/out.txt
./build/dev/usss/usss_driver usss/dataset/graph_ecommerce_updated GRR output/hnsn/out.txt
./build/dev/usss/usss_driver usss/dataset/graph_ecommerce_updated FLOW output/hnsn/out.txt
```

For the experiments above, the file `./output/hnsn/out.txt` contains one row per iteration for the MNP, IP (SuperGreedy++), and FW algorithms. Each row includes, in order: the iteration number, the density at that iteration, the norm of the load vector, and the time taken for that iteration (in nanoseconds, not cumulative). For the remaining seven methods, the file contains a single row each, representing iteration number 0, the final result produced by the algorithm, a placeholder value of 0, and the total time taken by the algorithm (in nanoseconds).

---

## Submodular Function Minimization (Contrapolymatroid Membership Problem)

### Overview 

Objective: Given a submodular function $f(S)$, we can define the base polymatroid $B(f)$ as a subspace of $R^{|V|}$. We would like to query whether arbitrary $y$ is in $B(f)$.

To do this, we notice that $y(S) - f(S) <= 0$, with $y(S)$ modular, so we can approximate an answer with SFM.

We implement the following algorithms:

- Frank-Wolfe 
- SuperGreedy++
- Minimum Norm Point (MNP)

## Directions to Generate Contrapolymatroid Dataset
We use the output of pushrelabel to generate our CM dataset. Specifically, `dsg/src/pushrelabel_contra.h` in DSG will write the final densest set to a file with the same name as the source dataset followed by a subscript `_contra`. This file is saved under the same directory as the source dataset.

Once all the datasets have been generated, go into `dss/src/perturbation.py` and point the file to the location of all datasets. Running it will create `<dataset_name>_perturbed.txt` files containing 5 contrapolymatroid query vectors, one on each line. The first corresponds to a YES instance, with the following ones being NO instances of decreasing difficulty corresponding to $\varepsilon = \{0.1, 1, 6, 12\}$ (see line 96 in `perturbation.py` to modify these).

For example, to generate a CM dataset for close-cliques, follow these steps:
1. Run the augmented push-relabel program (`pushrelabel_contra`) to extract the densest subgraph.
```bash
./build/dev/dsg/dsg_driver ./dsg/dataset/remap/close-cliques.txt pushrelabel_contra ./output/dsg/out.txt
```
Once the program terminates, a `./dsg/dataset/remap/close-cliques_contra.txt` file is created containing the densest subgraph.

2. Generate query vectors from the densest subgraph by running
```bash
python3 dss/src/perturbation.py
```
The above generates a `./dsg/dataset/remap/close-cliques_perturbed.txt` file containing the densest set, followed by a YES-instance vector, followed by 4 NO-instance vectors. Note that the script automatically looks for `_contra` suffixed files under `./dsg/dataset/remap/`.

Repeat step 1 for every base graph you'd like to create CM instances for. If several `_contra` suffixed files are contained under `remap`, a CM dataset will be created for each one of them in step 2. 

Our experiments reuse the DSS entrypoint but expect the index of the query vector (0-4) to be supplied. If no index is supplied, the 0 vector is used which recovers regular DSS.

### Running the Experiments

All algorithms are integrated via `dss/src/main.cpp`. While in the `root` directory:

```bash
./build/dev/dss/dss_driver <DATASET_PATH> <NAME_OF_ALGORITHM> <OUTPUT_PATH.txt> <ITERATIONS> <P_VALUE> <INDEX_OF_QUERY_VECTOR>
```

Note that the index of the query vector is $0$ for the YES instance, then $1,2,3,4$ for the NO instances corresponding to $\epsilon=1\mathrm{e}{-}1, 1, 6, 12$ respectively. 

### Examples 

```bash
./build/dev/dss/dss_driver ./dsg/dataset/remap/close-cliques.txt  frankwolfe ./output/cm/out.txt 10 1.1 2
./build/dev/dss/dss_driver ./dsg/dataset/remap/close-cliques.txt  mnp ./output/cm/out.txt 10 1.1 2
./build/dev/dss/dss_driver ./dsg/dataset/remap/close-cliques.txt  supergreedy ./output/cm/out.txt 10 1.1 2
```

The experiments above are using the 3rd query. Omitting the query will result in using the zero-vector (degenerates the regular p-mean DSG).

---

## Submodular Function Minimization (Minimum $s{-}t$ cut problem)

### Overview 

These are the experiments for the Minimum $s{-}t$ cut problem. We implement 4 algorithms:

- Frank-Wolfe 
- SuperGreedy++
- Minimum Norm Point (MNP)
- Edmonds Karp Flow-Based Algorithm for Minimum $s{-}t$ cut. 

### Datasets 
Available in the `sfm/dataset/` and `sfm/dataset/remap` directories:

- `sfm/dataset/m10.txt`
- `sfm/dataset/m30.txt`
- `sfm/dataset/m50.txt`
- `sfm/dataset/ww.txt`
- `sfm/src/dataset/BVZ-sawtooth/remap/BVZ-sawtooth0_s0_t1.max`
- `sfm/src/dataset/BVZ-sawtooth/remap/BVZ-sawtooth1_s0_t1.max`
- `sfm/src/dataset/BVZ-sawtooth/remap/BVZ-sawtooth2_s0_t1.max`
- `sfm/src/dataset/BVZ-sawtooth/remap/BVZ-sawtooth3_s0_t1.max`
- `sfm/src/dataset/BVZ-sawtooth/remap/BVZ-sawtooth4_s0_t1.max`
- `sfm/src/dataset/BVZ-sawtooth/remap/BVZ-sawtooth5_s0_t1.max`
- `sfm/src/dataset/BVZ-sawtooth/remap/BVZ-sawtooth6_s0_t1.max`
- `sfm/src/dataset/BVZ-sawtooth/remap/BVZ-sawtooth7_s0_t1.max`
- `sfm/src/dataset/BVZ-sawtooth/remap/BVZ-sawtooth8_s0_t1.max`
- `sfm/src/dataset/BVZ-sawtooth/remap/BVZ-sawtooth9_s0_t1.max`
- `sfm/src/dataset/BVZ-sawtooth/remap/BVZ-sawtooth10_s0_t1.max`
- `sfm/src/dataset/BVZ-sawtooth/remap/BVZ-sawtooth11_s0_t1.max`
- `sfm/src/dataset/BVZ-sawtooth/remap/BVZ-sawtooth12_s0_t1.max`
- `sfm/src/dataset/BVZ-sawtooth/remap/BVZ-sawtooth13_s0_t1.max`
- `sfm/src/dataset/BVZ-sawtooth/remap/BVZ-sawtooth14_s0_t1.max`
- `sfm/src/dataset/BVZ-sawtooth/remap/BVZ-sawtooth15_s0_t1.max`
- `sfm/src/dataset/BVZ-sawtooth/remap/BVZ-sawtooth16_s0_t1.max`
- `sfm/src/dataset/BVZ-sawtooth/remap/BVZ-sawtooth17_s0_t1.max`
- `sfm/src/dataset/BVZ-sawtooth/remap/BVZ-sawtooth18_s0_t1.max`
- `sfm/src/dataset/BVZ-sawtooth/remap/BVZ-sawtooth19_s0_t1.max`

### Running the Experiments

All algorithms are integrated via `sfm/src/main.cpp`. While in the `root` directory:

```bash
./build/dev/sfm/sfm_driver <DATASET_PATH> <NAME_OF_ALGORITHM> <OUTPUT_PATH.txt> <ITERATIONS_FOR_NONFLOW_ALGORITHMS>
```

### Examples 

```bash
./build/dev/sfm/sfm_driver ./sfm/dataset/m10.txt  frankwolfe ./output/mincut/out.txt 10
./build/dev/sfm/sfm_driver ./sfm/dataset/m10.txt  mnp ./output/cm/out.txt 10
./build/dev/sfm/sfm_driver ./sfm/dataset/m10.txt  supergreedy ./output/cm/out.txt 10
./build/dev/sfm/sfm_driver ./sfm/dataset/m10.txt  flow ./output/cm/out.txt 
```

Note that the value printed in the output file represents the negative of the min-cut; to obtain the actual min-cut value, you should negate it. Each row contains, in order: the iteration number, the negative min-cut value for that iteration, the norm of the load vector (for convex programming-based algorithms), and the time taken for that iteration (in nanoseconds, non-cumulative).

---

## Minimum Norm Point (MNP)

No additional datasets or experiments needed. Use the norm of the load vector statistics already collected across datasets and algorithms from DSG, DSS, USSS, and SFM tasks.