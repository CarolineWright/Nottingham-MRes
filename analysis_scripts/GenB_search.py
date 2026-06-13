#!/usr/bin/env python3

from Bio import Entrez, SeqIO
import os
import sys

# Set your email (NCBI requires this for API requests)
Entrez.email = "your.email@example.com"

# Specify the input file containing GenBank accession numbers
input_file = "Genbank_forTree.txt"

# Verify input file exists and is not empty
if not os.path.exists(input_file) or os.stat(input_file).st_size == 0:
    print(f"ERROR: Input file '{input_file}' is missing or empty!")
    sys.exit(1)

# Specify the output directory to save genome sequences
output_dir = "Genbank_Genomes"
os.makedirs(output_dir, exist_ok=True)

# Read accession numbers from the file
with open(input_file, "r") as f:
    accession_numbers = [line.strip() for line in f if line.strip()]

# Debugging: Check if we read accession numbers correctly
if not accession_numbers:
    print("ERROR: No accession numbers found in file.")
    sys.exit(1)

print(f"Downloading {len(accession_numbers)} genomes...")

# Loop through each accession number and download the genome
for acc in accession_numbers:
    try:
        print(f"Fetching {acc}...")
        handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
        record = handle.read()
        handle.close()

        # Verify content before writing to file
        if not record.strip():
            print(f"WARNING: No data retrieved for {acc}.")
            continue

        # Save to a FASTA file
        output_file = os.path.join(output_dir, f"{acc}.fasta")
        with open(output_file, "w") as f:
            f.write(record)
        
        print(f"Saved: {output_file}")
    
    except Exception as e:
        print(f"Error fetching {acc}: {e}")

# Verify output directory contents
if not os.listdir(output_dir):
    print("ERROR: No genomes were downloaded. Check accession numbers and network connection.")
else:
    print(f"Downloaded {len(os.listdir(output_dir))} genome files.")

print("All downloads completed.")
