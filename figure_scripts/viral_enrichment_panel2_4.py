#!/usr/bin/env python3
"""
Viral Enrichment Efficiency Plot
Compares viral read enrichment across different sequencing workflows
"""

import matplotlib.pyplot as plt
import numpy as np
import argparse


def create_viral_enrichment_plot(output_file='viral_enrichment_panel2_4.png'):

    pipelines = ['CL-SG', 'CL-CAP', 'CC-LYS', 'CC-PUR']
    pipeline_labels = ['Clinical\nShotgun', 'Clinical\nCapture', 'Culture\nLysate', 'Culture\nPurified']
    
    # Input data (% → fraction)
    # MODIFY THESE VALUES WITH YOUR DATA
    viral_percentages = [4.67E-04, 100.0, 3.01, 52.2]
    viral_fractions = [v / 100 for v in viral_percentages]

    # Reference (filtered = 1.0)
    filtered_fractions = [1.0] * 4

    color_schemes = {
        'CL-SG': plt.cm.Greys,
        'CL-CAP': plt.cm.Blues, 
        'CC-LYS': plt.cm.Reds,
        'CC-PUR': plt.cm.Greens
    }

    color_intensities = {
        'filtered': 0.5,
        'mapped': 0.7
    }

    bar_widths = {
        'filtered': 0.6,
        'mapped': 0.4
    }

    fig, ax = plt.subplots(figsize=(8, 6))
    x_pos = np.arange(len(pipelines))

    # Draw vertical nested bars
    for i, pipeline in enumerate(pipelines):
        cmap = color_schemes[pipeline]

        # Background (filtered = 1)
        ax.bar(i, filtered_fractions[i],
               width=bar_widths['filtered'],
               color=cmap(color_intensities['filtered']),
               alpha=0.8)

        # Foreground (viral fraction)
        ax.bar(i, viral_fractions[i],
               width=bar_widths['mapped'],
               color=cmap(color_intensities['mapped']),
               alpha=0.8)

    # Dynamic y-axis scale
    all_values = [v for v in viral_fractions if v > 0]
    min_val = min(all_values)

    ax.set_yscale('log')
    ax.set_ylim(min_val / 2, 1.05)

    # Labels
    ax.set_xticks(x_pos)
    ax.set_xticklabels(pipeline_labels)

    ax.set_ylabel('Fraction of filtered reads (log10 scale)',
                  fontsize=12, fontweight='bold')

    ax.set_title('Viral Enrichment Efficiency',
                 fontsize=14, fontweight='bold', pad=20)

    # Clean appearance
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, alpha=0.3, axis='y')

    ax.tick_params(axis='x', labelsize=10)
    ax.tick_params(axis='y', labelsize=10)

    # LOG-SAFE ANNOTATION
    LOG_OFFSET = 0.2  # constant shift in log space (~1.58×)
    REDUCED_LOG_OFFSET = LOG_OFFSET * 0.25  # 75% reduction
    
    y_max = 1.05

    for i, (pipeline, val) in enumerate(zip(pipelines, viral_fractions)):

        # Format label
        if val > 1e-3:
            label = f'{val:.3f}'
        else:
            label = f'{val:.2e}'

        # Apply reduced offset for selected pipelines
        if pipeline in ['CL-SG', 'CC-LYS', 'CC-PUR']:
            offset = REDUCED_LOG_OFFSET
        else:
            offset = LOG_OFFSET
        
        # Compute position (log-safe)
        y_text = val * (10 ** offset)
        

        # Prevent overflow
        if y_text > y_max * 0.9:
            y_text = val / (10 ** LOG_OFFSET)
            va = 'top'
        else:
            va = 'bottom'

        ax.text(i, y_text, label,
                ha='center', va=va,
                fontsize=9, fontweight='bold')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()

    # Summary (still in % for readability)
    print("\n=== Viral Enrichment Efficiency Results ===")
    print("Pipeline\t\t% Viral Reads\tFold vs CL-SG")
    print("-" * 50)

    baseline = viral_percentages[0]
    for pipeline, pct in zip(pipelines, viral_percentages):
        fold = pct / baseline if baseline > 0 else 0
        print(f"{pipeline}\t\t{pct:.2e}%\t{fold:.0f}×")

    return viral_fractions


def main():
    parser = argparse.ArgumentParser(description='Generate viral enrichment efficiency plot')
    parser.add_argument('--output', default='viral_enrichment_panel2_4.png')

    args = parser.parse_args()

    create_viral_enrichment_plot(args.output)

    print(f"Panel 2 figure saved as: {args.output}")


if __name__ == "__main__":
    main()
