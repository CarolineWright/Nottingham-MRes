#!/bin/bash

#SBATCH --job-name=bbmapDef_KX_job
#SBATCH --output=bbmapDef_KX_%j.out
#SBATCH --error=bbmapDef_KX_%j.err
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=480G

# Exit on error
set -e

# Define directories
output_dir_mapped="mappedKX"                  # Directory for mapped reads
output_dir_unmapped="unmappedKX"              # Directory for unmapped reads

# Reference file
REFERENCE="/path/to/reference/genome.fasta"

# Create output directories if they don't exist
mkdir -p "$output_dir_mapped"
mkdir -p "$output_dir_unmapped"

# Build the BBMap reference if not already done
bbmap.sh ref="$REFERENCE" nodisk

# Loop through all paired-end fastq.gz files in the current working directory
for r1_file in *_R1_filtered.fastq.gz; do
    # Extract base name without the R1/R2 suffix and extensions
    base_name=$(basename "$r1_file" "_R1_filtered.fastq.gz")
    r2_file="${base_name}_R2_filtered.fastq.gz"
    
    # Ensure the corresponding R2 file exists
    if [[ ! -f $r2_file ]]; then
        echo "Error: Paired R2 file for $r1_file not found. Skipping..."
        continue
    fi

    # Define output file paths
    mapped_output="$output_dir_mapped/${base_name}_mapped.fastq"
    unmapped_output="$output_dir_unmapped/${base_name}_unmapped.fastq"

    # Run BBMap for alignment
    echo "Aligning $r1_file and $r2_file..."
    bbmap.sh in="$r1_file" in2="$r2_file" outm="$mapped_output" outu="$unmapped_output" ref="$REFERENCE" nodisk
done

echo "Alignment completed. Mapped and unmapped reads are saved in their respective directories."
