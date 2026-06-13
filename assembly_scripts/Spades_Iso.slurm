#!/bin/bash

#SBATCH --job-name=SpIso
#SBATCH --output=SpIso_%j.out
#SBATCH --error=SpIso_%j.err
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=250G

INPUT_DIR="/path/to/interleaved/fastq/files"
OUTPUT_DIR="/path/to/output/spades"

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Loop through all interleaved .fastq.gz files in the input directory
for interleaved_file in "$INPUT_DIR"/*.fastq.gz; do
    # Get the base name of the file (without path and extension)
    base_name=$(basename "$interleaved_file" .fastq.gz)
    
    echo "Processing sample: $base_name"
    
    # Run SPAdes with the --isolate option and specify the output directory
    spades.py --12 "$interleaved_file" -o "$OUTPUT_DIR/$base_name" --isolate -k 21,33,55,77,99,127 --threads 10 --memory 250
    
    echo "Finished processing sample: $base_name"
done

echo "All samples processed."
