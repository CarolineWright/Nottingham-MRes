#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=8G
#SBATCH --job-name=FULL_RAxMLNG
#SBATCH --output=FULL_RAxMLNG.%j.out
#SBATCH --error=FULL_RAxMLNG.%j.err

# Set up Conda in a non-interactive environment
source /path/to/miniconda/etc/profile.d/conda.sh
conda deactivate 2>/dev/null
conda activate environment_core_analysis

###############################################################################
# Strict error handling
###############################################################################
set -euo pipefail

###############################################################################
# Ensure execution in the SLURM submission directory
###############################################################################
cd "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR not set}"

###############################################################################
# Basic logging and environment info
###############################################################################
echo "===== Job started at: $(date) ====="
echo "Hostname: $(hostname)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}"
echo "SLURM_JOB_CPUS_PER_NODE=${SLURM_JOB_CPUS_PER_NODE:-NA}"
echo "SLURM_MEM_PER_NODE=${SLURM_MEM_PER_NODE:-NA}"
echo "Working directory: $(pwd)"

###############################################################################
# Optional: record resource limits (useful for debugging OOM / scaling issues)
###############################################################################
ulimit -a

###############################################################################
# Define input and output files
###############################################################################
INPUT_ALIGNMENT="/path/to/alignment.fasta"
OUTPUT_PREFIX="analysis_prefix"

###############################################################################
# Run RAxML-NG: tree inference + bootstraps
###############################################################################
echo "===== Running RAxML-NG --all ====="
raxml-ng --all \
  --msa "$INPUT_ALIGNMENT" \
  --model GTR+G+F \
  --prefix "$OUTPUT_PREFIX" \
  --threads "${SLURM_CPUS_PER_TASK}"

###############################################################################
# Apply bootstrap support to the best tree
###############################################################################
echo "===== Applying bootstrap support to best tree ====="
raxml-ng --support \
  --tree "${OUTPUT_PREFIX}.raxml.bestTree" \
  --bs-trees "${OUTPUT_PREFIX}.raxml.bootstraps" \
  --prefix "${OUTPUT_PREFIX}_support"

###############################################################################
# Completion
###############################################################################
echo "===== Job finished successfully at: $(date) ====="
