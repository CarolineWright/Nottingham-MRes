#!/bin/bash

#SBATCH --job-name=fastp_filt
#SBATCH --output=fastp_filt_%j.out
#SBATCH --error=fastp_filt_%j.err
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=480G

# Input directory containing raw fastq.gz files
INPUT_DIR="/path/to/raw/fastq/files"

# Create the Filtered_fastq directory if it doesn't exist
mkdir -p Filtered_fastq2

# Loop through all paired-end fastq.gz files in the input directory
for file1 in "$INPUT_DIR"/*R1*.fastq.gz; do
    # Define corresponding R2 file
    file2="${file1/R1/R2}"
    
    # Check if R2 file exists
    if [[ -f "$file2" ]]; then
        # Set output file names
        base_name=$(basename "$file1" _R1.fastq.gz)  # Extract base name without _R1.fastq.gz
        output1="Filtered_fastq2/${base_name}_R1_filtered.fastq.gz"
        output2="Filtered_fastq2/${base_name}_R2_filtered.fastq.gz"
        json_output="Filtered_fastq2/${base_name}_fastp_stats.json"
        html_output="Filtered_fastq2/${base_name}_fastp_report.html"

        # Run fastp with additional trimming and filtering options
        fastp -i "$file1" -I "$file2" \
              -o "$output1" -O "$output2" \
              --dedup \
              --correction \
              --trim_poly_g \
              --length_required 110 \
              --json "$json_output" \
              --html "$html_output" \
              --thread 10

        echo "Processed: $file1 and $file2, outputs saved to $output1, $output2, $json_output, and $html_output"
    else
        echo "Warning: R2 file for $file1 not found. Skipping."
    fi
done

echo "fastp processing complete."
