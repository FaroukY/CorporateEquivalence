import os
import argparse

from config import load_config

"""
Generates slurms for the project in the specified config path.

Usage: python3 ./configs/generate_slurms --config <config_path>
Ex: python3 ./configs/generate_slurms.py --config ./configs/dsg.yaml
"""


def generate_slurm_script(
    config, entrypoint, output_path, args, algorithm_name, dataset_path
):
    args_str = " ".join(map(str, args))
    script = f"""#!/bin/bash
#SBATCH --mem={config.mem} # Memory allocation
#SBATCH --nodes={config.nodes} # Number of nodes
#SBATCH --cpus-per-task={config.cpus_per_task} # Number of CPU cores per task
#SBATCH --tasks={config.tasks}
#SBATCH --partition={config.partition} # Partition name
#SBATCH --account={config.account} # Account/project name
#SBATCH --job-name={config.job_name} # Job name
#SBATCH --time={config.time} # Time limit (HH:MM:SS)
#SBATCH --error={config.error} # Standard error log
#SBATCH --output={config.output} # Standard output log

# ---------------------------
# Run the C++ Script (from root)
# ---------------------------
module load intel-tbb/2021.9.0
module load gurobi/10.0.1
export GRB_LICENSE_FILE=/sw/external/gurobi/gurobi.lic
{entrypoint} {dataset_path} {algorithm_name} {output_path} {args_str}
"""
    return script


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate experiment slurms")
    parser.add_argument("--config", type=str)
    parser.add_argument("--slurm_dir", type=str, default="slurms")
    parser.add_argument("--project", type=str, default="")

    args_cmd = parser.parse_args()
    args = load_config(
        config_path=os.path.abspath(args_cmd.config),
    )

    # If no project is specified, use the config name
    args.project = (
        os.path.splitext(os.path.basename(args_cmd.config))[0]
        if args_cmd.project == ""
        else args_cmd.project
    )
    args.slurm_dir = args_cmd.slurm_dir

    # Create the save directory if it doesn't exist
    save_dir = f"{args.slurm_dir}/{args.project}"
    os.makedirs(save_dir, exist_ok=True)

    # Create the slurm scripts
    for dataset, dataset_config in args.datasets.items():
        for algorithm, experiments in args.algorithms.items():
            for idx, algorithm_config in enumerate(experiments.experiments):
                experiment_name = algorithm_config.get("name", idx)
                config = args.default.copy()
                # Merge in dataset overrides, followed by algorithm overrides (which take precedence)
                if dataset_config.get("overrides") is not None:
                    for k, v in dataset_config.overrides.items():
                        setattr(config, k, v)
                if algorithm_config.get("overrides") is not None:
                    for k, v in algorithm_config["overrides"].items():
                        setattr(config, k, v)
                # Generate the slurm script
                slurm_script = generate_slurm_script(
                    config=config,
                    entrypoint=args.entrypoint,
                    output_path=f"{args.output_dir}/{dataset}_{algorithm}_{experiment_name}.txt",
                    args=algorithm_config["args"],
                    algorithm_name=algorithm,
                    dataset_path=dataset_config.dataset_path,
                )

                # Write the slurm script to a file
                slurm_filename = f"{args.slurm_dir}/{args.project}/{dataset}_{algorithm}_{experiment_name}.slurm"
                with open(slurm_filename, "w") as f:
                    f.write(slurm_script)
                print(f"Generated slurm script: {slurm_filename}")
    print(f"All slurm scripts generated in {save_dir}")
