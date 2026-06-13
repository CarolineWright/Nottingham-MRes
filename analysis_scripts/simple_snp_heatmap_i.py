#!/usr/bin/env python3
"""
Simple SNP Distance Heatmap Generator
Creates individual heatmap for a single FASTA alignment
Counts gaps as genuine differences (important for indel analysis)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import sys

def calculate_snp_distances(fasta_file):
    """
    Calculate pairwise SNP distances from alignment
    Counts gaps as differences (for indel analysis)
    """
    print(f"Reading sequences from {fasta_file}...")
    
    # Read sequences
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq).upper()
    
    seq_names = list(sequences.keys())
    print(f"Found {len(seq_names)} sequences:")
    for name in seq_names:
        print(f"  - {name}")
    
    # Calculate distances
    n_seqs = len(seq_names)
    distance_matrix = np.zeros((n_seqs, n_seqs), dtype=int)
    
    print("Calculating pairwise distances...")
    for i in range(n_seqs):
        for j in range(i, n_seqs):
            if i == j:
                distance_matrix[i][j] = 0
            else:
                seq1 = sequences[seq_names[i]]
                seq2 = sequences[seq_names[j]]
                
                # Count differences (including gaps vs bases)
                differences = 0
                valid_sites = 0
                
                for pos in range(min(len(seq1), len(seq2))):
                    base1, base2 = seq1[pos], seq2[pos]
                    
                    # Skip position if any base in ignore_chars set
                    ignore_chars = set('?')   # previously: set('N?')
                    if base1 not in ignore_chars and base2 not in ignore_chars:
                        valid_sites += 1
                        if base1 != base2:  
                            differences += 1
                
                distance_matrix[i][j] = differences
                distance_matrix[j][i] = differences
    
    # Convert to DataFrame
    distance_df = pd.DataFrame(distance_matrix, 
                              index=seq_names, 
                              columns=seq_names)
    
    return distance_df

def clean_sample_names(names):
    """Clean up sequence names for better display"""
    cleaned = []
    for name in names:
        # Remove common prefixes and simplify names
        clean = name.replace('LSDV_', '').replace('Mongolia_', 'M')
        
        # Shorten years
        clean = clean.replace('_2022', '_22').replace('_2021', '_21')
        
        # Highlight assembly differences
        if '_Haga24' in clean:
            clean = clean.replace('_Haga24', '_H')
        elif '_i' in clean:
            clean = clean.replace('_i', '_I')
            
        cleaned.append(clean)
    
    return cleaned

def create_heatmap(distance_df, output_file, title="Pairwise SNP Distances"):
    """Create publication-quality heatmap"""
    
    # Clean labels
    clean_labels = clean_sample_names(distance_df.index.tolist())
    
    # Set up plot
    plt.figure(figsize=(10, 8))
    
    # Create heatmap with custom styling
    ax = sns.heatmap(distance_df, 
                     annot=True,           # Show values
                     fmt='d',              # Integer format
                     cmap='RdYlBu_r',      # Red-Yellow-Blue reversed
                     square=True,          # Square cells
                     linewidths=0.5,       # Cell borders
                     cbar_kws={'label': 'SNP Distance', 'shrink': 0.8},
                     xticklabels=clean_labels,
                     yticklabels=clean_labels,
                     annot_kws={'size': 11})  # Annotation font size
    
    # Styling
    plt.title(title, fontsize=16, fontweight='bold', pad=20)
    plt.xlabel('Isolates', fontsize=12, fontweight='bold')
    plt.ylabel('Isolates', fontsize=12, fontweight='bold')
    
    # Rotate labels
    plt.xticks(rotation=45, ha='right', fontsize=10)
    plt.yticks(rotation=0, fontsize=10)
    
    # Save with high quality
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    
    return distance_df

def print_summary_stats(distance_df, dataset_name):
    """Print summary statistics for the distance matrix"""
    print(f"\n=== {dataset_name} Summary Statistics ===")
    
    # Get upper triangle (excluding diagonal) for unique pairwise distances
    upper_triangle = distance_df.values[np.triu_indices_from(distance_df.values, k=1)]
    
    print(f"Number of sequences: {len(distance_df)}")
    print(f"Number of pairwise comparisons: {len(upper_triangle)}")
    print(f"Mean distance: {upper_triangle.mean():.2f} SNPs")
    print(f"Maximum distance: {upper_triangle.max()} SNPs")
    print(f"Minimum distance: {upper_triangle.min()} SNPs")
    
    if len(upper_triangle) > 0:
        print(f"Distance range: {upper_triangle.min()} - {upper_triangle.max()} SNPs")
    
    # Show the actual distance matrix
    print(f"\nDistance Matrix:")
    print(distance_df.to_string())

def main():
    """Main function"""
    if len(sys.argv) != 3:
        print("Usage: python simple_snp_heatmap_i.py <fasta_file> <output_prefix>")
        print("Example: python simple_snp_heatmap_i.py alignment.fasta figure_output")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    output_prefix = sys.argv[2]
    
    # Determine dataset type for title
    if "Mine" in fasta_file:
        title = "Isolates - Your Assemblies"
        dataset_name = "Your Assemblies"
    elif "Haga" in fasta_file:
        title = "Isolates - Reference Assemblies"  
        dataset_name = "Reference Assemblies"
    elif "ALL" in fasta_file or "Combined" in fasta_file:
        title = "Isolates - Combined Assemblies"
        dataset_name = "Combined Assemblies"
    else:
        title = "Isolates - SNP Distances"
        dataset_name = "Unknown Dataset"
    
    # Calculate distances
    distance_df = calculate_snp_distances(fasta_file)
    
    # Create heatmap
    output_png = f"{output_prefix}.png"
    create_heatmap(distance_df, output_png, title)
    print(f"Heatmap saved: {output_png}")
    
    # Save distance matrix as CSV
    output_csv = f"{output_prefix}_distances.csv"
    distance_df.to_csv(output_csv)
    print(f"Distance matrix saved: {output_csv}")
    
    # Print summary statistics
    print_summary_stats(distance_df, dataset_name)

if __name__ == "__main__":
    main()
