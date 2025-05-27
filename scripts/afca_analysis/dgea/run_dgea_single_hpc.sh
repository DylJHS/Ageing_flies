#!/bin/bash

# This SLURM batch script submits an array job to perform differential gene expression analysis (DGEA) across multiple cell types in parallel. 
# Each array task processes a single .h5ad file using run_dgea_single.R, which is run via Rscript. 
# Input files are read from a specified directory, and results are saved to a corresponding output folder. 
# The job requests 4 CPUs, 16 GB of memory, and includes email notifications and logging for each task.

#SBATCH --job-name=dgea_single
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/ageing_flies/logs/dgea_%A_%a.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/ageing_flies/logs/dgea_%A_%a.err
#SBATCH --time=05:00:00
#SBATCH --array=
#SBATCH --mem=32G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

set -euo pipefail  

# Load necessary modules and activate conda environment
export PATH=/hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/bin:$PATH
source /hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/etc/profile.d/conda.sh
conda activate /hpc/shared/onco_janssen/dhaynessimmons/envs/dgea_env

# Paths
INPUT_DIR="/hpc/shared/onco_janssen/dhaynessimmons/projects/ageing_flies/data/ct_specific"
OUTPUT_DIR="/hpc/shared/onco_janssen/dhaynessimmons/projects/ageing_flies/results/afca_analysis/dgea_results/ct_specific"

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Build array of files
FILES=($INPUT_DIR/*)
INPUT_FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

# Run the per-file R script
Rscript /hpc/shared/onco_janssen/dhaynessimmons/projects/ageing_flies/scripts/afca_analysis/dgea/run_dgea_single_hpc.r "$INPUT_FILE" "$OUTPUT_DIR"
