#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --job-name=3seq_COR
#SBATCH --output=3seq_COR.%j.out
#SBATCH --error=3seq_COR.%j.err

###############################################################################
# SLURM Script: Run 3seq recombination detection analysis
###############################################################################

###############################################################################
# Initialize environment
###############################################################################
eval "$(micromamba shell hook --shell bash)"
micromamba deactivate 2>/dev/null || true
conda deactivate 2>/dev/null || true
micromamba activate environment_recombination

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
echo "===== 3seq Recombination Analysis Started at: $(date) ====="
echo "Hostname: $(hostname)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}"
echo "SLURM_MEM_PER_NODE=${SLURM_MEM_PER_NODE:-NA}"
echo "Working directory: $(pwd)"
echo ""

###############################################################################
# Display resource limits
###############################################################################
echo "===== Resource Limits ====="
ulimit -a
echo ""

###############################################################################
# Define input and output files
###############################################################################
INPUT_ALIGNMENT="/path/to/alignment.fasta"
OUTPUT_PREFIX="analysis_prefix"

###############################################################################
# Run 3seq: Recombination Detection
###############################################################################
echo "===== Running 3seq fullrun mode ====="
echo "Input file: $INPUT_ALIGNMENT"
echo "Output prefix: $OUTPUT_PREFIX"
echo "Options: -d (distinct sequences only), -t (p-value threshold)"
echo ""

# Pipe 'y' to automatically answer confirmation prompts
echo "y" | 3seq -f "$INPUT_ALIGNMENT" \
     -d \
     -id "$OUTPUT_PREFIX" \
     -t0.05

###############################################################################
# Verify output files were created
###############################################################################
echo ""
echo "===== Verifying Output Files ====="
if [ -f "${OUTPUT_PREFIX}.3s.rec" ]; then
    echo "✓ Main results file created: ${OUTPUT_PREFIX}.3s.rec"
    REC_LINES=$(wc -l < "${OUTPUT_PREFIX}.3s.rec")
    echo "  Number of recombinant triplets: $((REC_LINES - 1))"
else
    echo "Results file not found (may indicate no recombination detected)"
fi

if [ -f "${OUTPUT_PREFIX}.3s.log" ]; then
    echo "✓ Log file created: ${OUTPUT_PREFIX}.3s.log"
fi

echo ""
echo "===== 3seq Recombination Analysis Completed at: $(date) ====="
