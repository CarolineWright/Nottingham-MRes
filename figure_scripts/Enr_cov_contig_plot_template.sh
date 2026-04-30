#!/bin/bash
#SBATCH --job-name=lsdv_enr_plot
#SBATCH --output=lsdv_enr_plot_%j.out
#SBATCH --error=lsdv_enr_plot_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

#=============================================================================
# LSDV Coverage + Contig Visualization Script
# Author: Caroline Wright
# Purpose: Generate stacked coverage and contig alignment plots
# FILE: Nottingham-MRes/figure_scripts/Enr_cov_contig_plot_template.sh
#=============================================================================


#=============================================================================
# CONFIGURATION - Edit these variables for your system
#=============================================================================
BASE_DIR="/path/to/your/project"
COVERAGE_DIR="${BASE_DIR}/Coverage_plots/Combined_regions"  
CONTIG_DIR="${BASE_DIR}/Spades_assemblies"
OUTPUT_DIR="${BASE_DIR}/figures"
CONDA_PATH="/path/to/conda"
CONDA_ENV="biofigures"

# File specifications
COVERAGE_FILE="${COVERAGE_DIR}/Enr_regions.tsv"
CONTIG_FILE="${CONTIG_DIR}/Enr_contig_alignments.tsv"  
OUTPUT_PLOT="${OUTPUT_DIR}/Enr_cov_contig_plot.png"


# Setup conda
source "${CONDA_PATH}/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"


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
fig, ax = plt.subplots(figsize=(16,6))  # Increased height for more rows

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

    # Adjust x positions to center 100bp bins on their midpoints
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

    coverage_offset += 1.2  # Space between coverage tracks

###############################################################################
# Contig tracks (below axis) - Improved multi-row layout with jitter
###############################################################################
line_segments = []
segment_colors = []

current_y = -1

for sample in samples:
    df = contigs[contigs["sample"] == sample]
    
    if df.empty:
        print(f"Warning: no contig data found for sample {sample}")
        continue
    
    num_contigs = len(df)
    print(f"Sample {sample}: {num_contigs} contigs")
    
    # Calculate rows needed (aim for ~15-20 contigs per row max)
    if num_contigs > 100:
        rows_for_sample = min(8, max(4, num_contigs // 15))  # 4-8 rows for highly fragmented
    elif num_contigs > 30:
        rows_for_sample = min(4, max(2, num_contigs // 15))  # 2-4 rows for moderately fragmented  
    else:
        rows_for_sample = 1  # Single row for low fragmentation
    
    print(f"Sample {sample}: Using {rows_for_sample} rows")
    
    # Sample-specific row spacing
    if sample == "05":  # Sample 05 - tighter spacing for many contigs
        row_spacing = 0.25
    elif sample == "07":  # Sample 07 - wider spacing to avoid overlap of longest contig
        row_spacing = 0.35
    else:
        row_spacing = 0.3  # Default spacing


    # Add sample label (centered vertically for multi-row samples)
    sample_height = rows_for_sample * 0.2 + (rows_for_sample-1) * 0.1  # Row height + spacing
    label_y_position = current_y - sample_height/2
    ax.text(
        -2000, 
        label_y_position, 
        sample, 
        ha='right', 
        va='center', 
        fontsize=10
    )
    
    # Sort contigs by start position for better visualization
    df_sorted = df.sort_values('start')
    starts = df_sorted["start"].to_numpy()
    ends = df_sorted["end"].to_numpy()
    
    # Use default color if sample not in colors dict
    sample_color = colors.get(sample, "#228B22")

    # Distribute contigs across rows with some randomization to avoid perfect stacking
    for idx, (s, e) in enumerate(zip(starts, ends)):
        row_num = idx % rows_for_sample
        row_offset = row_num * row_spacing  # Use sample-specific spacing
        
        # Add small random jitter within each row
        jitter = np.random.uniform(-0.08, 0.08)
        y_level = current_y - row_offset + jitter
        
        line_segments.append([(s, y_level), (e, y_level)])
        segment_colors.append(sample_color)
    
    # Move to next sample position
    current_y -= (rows_for_sample * row_spacing + 0.5)  # Space between samples

# CREATE AND ADD THE LINE COLLECTION (this is what you deleted)
lc = LineCollection(
    line_segments,
    colors=segment_colors,
    linewidths=1.5  # No transparency - removed alpha parameter
)

ax.add_collection(lc)

###############################################################################
# Adjust vertical limits to ensure contigs remain visible
###############################################################################
ax.set_ylim(current_y - 0.5, coverage_offset + 2)

###############################################################################
# Add ITR labels
###############################################################################
ax.text(
    (left_itr_start + left_itr_end) / 2,
    coverage_offset + 0.2,
    "ITR",
    ha="center",
    fontsize=9
)

ax.text(
    (right_itr_start + right_itr_end) / 2,
    coverage_offset + 0.2,
    "ITR",
    ha="center",
    fontsize=9
)

###############################################################################
# Axis styling
###############################################################################
ax.set_xlim(1, genome_length)

ax.set_xlabel("Genome position")
ax.set_yticks([])

ax.spines["left"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_visible(False)

###############################################################################
# Reference label
###############################################################################
ax.text(
    genome_length/2,
    -0.15,
    "LSDV_Israel_2012_KX894508",
    ha="center",
    va="top",
    fontsize=11
)

###############################################################################
# Save plot
###############################################################################
plt.tight_layout()
plt.savefig(output_plot, dpi=300, bbox_inches="tight")

print("Plot saved to:", output_plot)

EOF
