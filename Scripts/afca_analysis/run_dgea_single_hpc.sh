# run_dgea_array.sh
#!/bin/bash
#SBATCH --job-name=dgea_single
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/logs/dgea_%A_%a.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/logs/dgea_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=0-58%5
#SBATCH --mem=16G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

module load R

# Paths
INPUT_DIR="/hpc/shared/onco_janssen/dhaynessimmons/data/ageing_flies/afca_data/base_body_h5ad_data/ct_specific"
OUTPUT_DIR="/hpc/shared/onco_janssen/dhaynessimmons/results/ageing_flies/afca_analysis/dgea_results/ct_specific"

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Build array of files
FILES=($INPUT_DIR/*)
INPUT_FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

# Run the per-file R script
Rscript run_dgea_single.R "$INPUT_FILE" "$OUTPUT_DIR"
