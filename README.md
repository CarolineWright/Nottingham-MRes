# Nottingham-MRes

**University:** University of Nottingham  
**School:** School of Veterinary Medicine and Science  
**Course Title:** Level 7 Bioinformatics Scientist (MRes) Degree Apprenticeship  
**Organisation/Employer:** The Pirbright Institute  

## Description

This repository contains supplementary data and scripts (to follow) associated with my Master's research project. It includes:

- `supplementary_data/`  
  - Combined PDF (Supplementary Data 1 & 2)  
  - Tabulated data in CSV format (Supplementary Data 3 & 4)  
  - *Note:* Supplementary Data 3 was filtered using sequence identity (≥70%) and query coverage (≥70%) thresholds, retaining only the top hit per query (lowest E-value).

- `qc_scripts/`  
  Scripts for preprocessing and quality control of raw sequencing reads.

- `assembly_scripts/`  
  Scripts for genome assembly and polishing.

- `analysis_scripts/`  
  Scripts for consensus generation, variant calling, alignment, BLASTp alignmnet and result filtering/tabulation etc.

- `figure_scripts/`  
  Scripts for figure and plot generation.

## Requirements

- Python 3.10  
*(Dependencies for individual scripts will be specified within each script or in an accompanying requirements file.)*

## Acknowledgements

I would like to thank Dr Tim Downing (The Pirbright Institute) for their contribution to this MRes thesis project, particularly in the area of recombination analysis and associated scripting.  
Thanks also to the members of the High Throughput Sequencing Group (The Pirbright Institute) for operating the Illumina sequencing facility that generated the data presented in this thesis.  

I would also like to acknowledge the funders and fellow project participants of BBSRC project BB/Y007409/1, titled *"Understanding and controlling lumpy skin disease in Indian cattle"*, who supported the work associated with Supplementary Data 3 & 4.
