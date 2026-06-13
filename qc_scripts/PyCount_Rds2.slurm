#!/bin/bash

#SBATCH --job-name=read_cnt
#SBATCH --output=read_cnt_%j.out
#SBATCH --error=read_cnt_%j.err
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=480G

# Exit on error
set -e

# Activate appropriate virtual environment
source /path/to/miniconda/etc/profile.d/conda.sh && conda activate environment_python_analysis

# Create a Python script for counting reads
cat << EOF > read_cnt_script.py
import os
import gzip
import re

# Function to count the number of reads in a FASTQ file
def count_reads_in_fastq(file_path):
    with gzip.open(file_path, 'rt') as f:
        count = sum(1 for _ in f) // 4  # Each read spans 4 lines
    return count

# Function to process paired-end FASTQ files and calculate total reads
def process_fastq_files(output_file="read_counts.txt"):
    fastq_files = [f for f in os.listdir('.') if f.endswith('.fastq.gz')]
    
    # Group files by base name (assuming _R1 and _R2 naming convention)
    paired_files = {}
    for file in fastq_files:
        match = re.match(r"(.+)_R[12]_filtered.fastq.gz", file)
        if match:
            base_name = match.group(1)
            paired_files.setdefault(base_name, []).append(file)

    with open(output_file, 'w') as out:
        out.write("Base Name\tTotal Reads\n")
        for base_name, files in paired_files.items():
            if len(files) == 2:  # Ensure both R1 and R2 exist
                r1_file = next(f for f in files if "_R1" in f)
                r2_file = next(f for f in files if "_R2" in f)

                # Count reads in both files
                r1_count = count_reads_in_fastq(r1_file)
                r2_count = count_reads_in_fastq(r2_file)

                # Verify counts match (optional sanity check)
                if r1_count != r2_count:
                    print(f"Warning: Read counts differ for {base_name}: R1={r1_count}, R2={r2_count}")

                total_reads = r1_count + r2_count
                out.write(f"{base_name}\t{total_reads}\n")

            else:
                print(f"Warning: Missing paired-end files for {base_name}")

if __name__ == "__main__":
    process_fastq_files()
EOF

# Run the Python script
python3 read_cnt_script.py

# Clean up the Python script
rm -f read_cnt_script.py
