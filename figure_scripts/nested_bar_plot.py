#!/usr/bin/env python3
"""
Sequencing Efficiency Nested Bar Plot Generator
Compares read retention across different library preparation strategies
"""

import matplotlib.pyplot as plt
import numpy as np
import argparse

def parse_sequencing_data(input_file):
    """Parse tab-delimited input file with sequencing metrics"""
    data = {}
    
    with open(input_file, 'r') as f:
        header = f.readline().strip().split('\t')
        
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                pipeline = parts[0]
                raw = int(parts[1])
                filtered = int(parts[2])
                mapped = int(parts[3])
                unique = int(parts[4])
                
                data[pipeline] = {
                    'raw': raw,
                    'filtered': filtered,
                    'mapped': mapped,
                    'unique': unique
                }
    
    return data

def create_nested_bar_plot(data, output_file='sequencing_efficiency.png'):
    """Create nested bar plot showing read retention across pipelines"""
    
    pipeline_order = ['CL-SG', 'CL-CAP', 'CC-LYS', 'CC-PUR']
    pipeline_labels = ['Clinical\nShotgun', 'Clinical\nCapture', 'Culture\nLysate', 'Culture\nPurified']
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    y_pos = np.arange(len(pipeline_order))
    
    color_schemes = {
        'CL-SG': plt.cm.Greys,
        'CL-CAP': plt.cm.Blues,
        'CC-LYS': plt.cm.Reds,
        'CC-PUR': plt.cm.Greens
    }
    
    color_intensities = {
        'raw': 0.4,
        'filtered': 0.6,
        'mapped': 0.8,
        'unique': 0.9
    }
    
    proportions = {}
    
    for i, pipeline in enumerate(pipeline_order):
        if pipeline in data:
            raw_reads = data[pipeline]['raw']
            filtered_reads = data[pipeline]['filtered']
            mapped_reads = data[pipeline]['mapped']
            unique_reads = data[pipeline]['unique']
            
            proportions[pipeline] = {
                'raw': 100.0,
                'filtered': (filtered_reads / raw_reads * 100) if raw_reads > 0 else 0,
                'mapped': (mapped_reads / raw_reads * 100) if raw_reads > 0 else 0,
                'unique': (unique_reads / raw_reads * 100) if raw_reads > 0 else 0
            }
    
    bar_widths = {
        'raw': 0.8,
        'filtered': 0.7,
        'mapped': 0.6,
        'unique': 0.9
    }
    
    # Create nested bars for each pipeline
    for i, pipeline in enumerate(pipeline_order):
        if pipeline in proportions:
            cmap = color_schemes[pipeline]
            
            # Draw bars from largest (background) to smallest (foreground)
            stages = ['raw', 'filtered', 'mapped', 'unique']
            
            for stage in stages:
                width = proportions[pipeline][stage]
                color = cmap(color_intensities[stage])
                bar_width = bar_widths[stage]

                ax.barh(i, width, height=bar_width, color=color, alpha=0.8)
    
    # Customize the plot
    ax.set_yticks(y_pos)
    ax.set_yticklabels(pipeline_labels)
    ax.tick_params(axis='y', labelsize=12)
    ax.set_xlabel('Retention as Proportion of Raw Reads (%)', fontsize=16, fontweight='bold')
    ax.set_title('Sequencing Read Retention Across Library Preparation Strategies', 
                 fontsize=14, fontweight='bold', pad=20)
    
    # Remove top/right spines for clean look
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Set x-axis limits and scale
    ax.set_xscale('log')
    ax.set_xlim(0.001, 105)
    
    # Add grid for readability
    ax.grid(True, alpha=0.3, axis='x')
    
    # Add percentage labels for unique reads (most important metric)
    for i, pipeline in enumerate(pipeline_order):
        if pipeline in proportions:
            unique_pct = proportions[pipeline]['unique']
            if unique_pct > 0.01:  # Only label if visible
                ax.text(unique_pct + 0.1, i, f'{unique_pct:.2f}%', 
                       va='center', ha='left', fontsize=9, fontweight='bold')
            else:
                ax.text(0.002, i, f'{unique_pct:.4f}%', 
                       va='center', ha='left', fontsize=9, color='red')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    
    # Print summary statistics
    print("\n=== Sequencing Efficiency Summary ===")
    print("Pipeline\t\tFiltering\tViral Enrich\tUnique Yield")
    print("-" * 60)
    
    for pipeline in pipeline_order:
        if pipeline in proportions:
            filtering_eff = proportions[pipeline]['filtered']
            viral_enrich = proportions[pipeline]['mapped'] / proportions[pipeline]['filtered'] * 100 if proportions[pipeline]['filtered'] > 0 else 0
            unique_yield = proportions[pipeline]['unique']
            
            print(f"{pipeline}\t\t{filtering_eff:.2f}%\t\t{viral_enrich:.3f}%\t\t{unique_yield:.4f}%")
    
    return proportions

def main():
    parser = argparse.ArgumentParser(description='Generate nested bar plot for sequencing pipeline efficiency')
    parser.add_argument('input_file', help='Tab-delimited input file with sequencing data')
    parser.add_argument('--output', default='sequencing_efficiency.png', help='Output PNG file')
    
    args = parser.parse_args()
    
    print(f"Reading data from {args.input_file}...")
    data = parse_sequencing_data(args.input_file)
    
    print("Creating nested bar plot...")
    proportions = create_nested_bar_plot(data, args.output)
    
    print(f"Plot saved as {args.output}")

if __name__ == "__main__":
    main()
