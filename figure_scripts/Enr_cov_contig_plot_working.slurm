#!/bin/bash
#SBATCH --job-name=lsdv_enr_plot
#SBATCH --output=lsdv_enr_plot_%j.out
#SBATCH --error=lsdv_enr_plot_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

# Set up Conda in a non-interactive environment
source /path/to/miniconda/etc/profile.d/conda.sh
conda deactivate 2>/dev/null
conda activate environment_visualization

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
# USER INPUT - MODIFY THESE PATHS
###############################################################################
COVERAGE_DIR="/path/to/coverage/data"
CONTIG_DIR="/path/to/contig/alignments"
OUTPUT_DIR="/path/to/output/figures"

mkdir -p "$OUTPUT_DIR"

COVERAGE_FILE="$COVERAGE_DIR/coverage_data.tsv"
CONTIG_FILE="$CONTIG_DIR/contig_alignments.tsv"

OUTPUT_PLOT="$OUTPUT_DIR/coverage_contig_plot.pdf"

###############################################################################
# Run Python plotting
###############################################################################
python << EOF

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

###############################################################################
# Input files
###############################################################################
coverage_file = "$COVERAGE_FILE"
contig_file = "$CONTIG_FILE"
output_plot = "$OUTPUT_PLOT"

###############################################################################
# Genome parameters
###############################################################################
genome_length = 150773

left_itr_start = 1
left_itr_end = 2418

right_itr_start = 148000
right_itr_end = 150773

###############################################################################
# Load data
###############################################################################
cov = pd.read_csv(coverage_file, sep="\t")
contigs = pd.read_csv(contig_file, sep="\t")

# Keep sample names as strings and ensure consistent formatting
cov["sample"] = cov["sample"].astype(str).str.zfill(2)  # This will pad to 2 digits
contigs["sample"] = contigs["sample"].astype(str).str.zfill(2)

samples = sorted(cov["sample"].unique())
print("Samples detected:", samples)

cov = cov[cov["sample"].isin(samples)]
contigs = contigs[contigs["sample"].isin(samples)]

cov.sort_values(["sample","start"], inplace=True)
contigs.sort_values(["sample","start"], inplace=True)

###############################################################################
# Colormap colours
###############################################################################
cmap = plt.cm.Blues

colors = {
    "05": cmap(0.4),
    "07": cmap(0.7)
}

###############################################################################
# Create figure
###############################################################################
fig, ax = plt.subplots(figsize=(16,6))

###############################################################################
# Draw solid reference axis
###############################################################################
ax.hlines(
    y=0,
    xmin=1,
    xmax=genome_length,
    linewidth=1,
    color="black"
)

###############################################################################
# Shade ITR regions
###############################################################################
ax.axvspan(
    left_itr_start,
    left_itr_end,
    color="grey",
    alpha=0.2,
    zorder=0
)

ax.axvspan(
    right_itr_start,
    right_itr_end,
    color="grey",
    alpha=0.2,
    zorder=0
)

###############################################################################
# Coverage tracks (above axis)
###############################################################################
coverage_offset = 0

for sample in samples:

    df = cov[cov["sample"] == sample]

    if df.empty:
        print(f"Warning: no coverage data found for sample {sample}")
        continue

    x = df["start"].to_numpy() + 50  # Center each bin on its midpoint
    y = df["coverage"].to_numpy()

    if len(y) == 0:
        print(f"Warning: coverage array empty for sample {sample}")
        continue

    # normalize coverage
    y = y / y.max()

    ax.fill_between(
        x,
        coverage_offset,
        y + coverage_offset,
        color=colors[sample],
        linewidth=0.5,
        edgecolor=colors[sample]
    )

    # Add sample label next to coverage track
    ax.text(
        -2000,
        coverage_offset + 0.5,
        sample,
        ha='right',
        va='center',
        fontsize=10
    )

    coverage_offset += 1.2

###############################################################################
# Contig visualization (below axis)
###############################################################################

current_y = -1

for sample in samples:
    df = contigs[contigs["sample"] == sample]
    
    if df.empty:
        print(f"Warning: no contig data found for sample {sample}")
        continue
    
    num_contigs = len(df)
    print(f"Sample {sample}: {num_contigs} contigs")

    ax.text(
        -2000, 
        current_y, 
        sample, 
        ha='right', 
        va='center', 
        fontsize=10
    )
    
    for idx, row in df.iterrows():
        x_start = row['start']
        x_end = row['end']
        length = row['length']
        
        ax.hlines(
            y=current_y,
            xmin=x_start,
            xmax=x_end,
            color=colors[sample],
            linewidth=3
        )
        
        current_y -= 0.8

###############################################################################
# Formatting
###############################################################################
ax.set_xlim(-3000, genome_length + 3000)
ax.set_ylabel('Coverage & Contigs', fontsize=12)
ax.set_xlabel('Genome Position (bp)', fontsize=12)
ax.set_title('Coverage and Contig Alignment', fontsize=14, fontweight='bold')

plt.tight_layout()
plt.savefig(output_plot, dpi=300, bbox_inches='tight')
plt.show()

print(f"Plot saved to: {output_plot}")

EOF

echo "===== Job finished successfully at: $(date) ====="
