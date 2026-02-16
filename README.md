# Description

This repository contains the analysis code used in the manuscript entitled "Mechano-chemosensitive regulation of T cell migration by Filamin A via CXCR4 and β1 integrin dynamic remodeling", submitted to Nature Communications.

The repository is organized into three independent analysis pipelines corresponding to the different experimental approaches described in the study.

# Analysis Pipelines

## 1. T cell contact area quantification (folder: T_cell_contact_area)
- Image processing and quantification of T cell contact area measurements.
- Scripts include segmentation and quantification of binary images. 
- Quantitative measurements were exported as CSV files for downstream statistical analysis. 
- Graphs and statistical analyses were generated using GraphPad Prism.

## 2. Single-molecule TIRF analysis (folder: Single_molecule_TIRF)
- Analysis of CXCR4 and β1 integrin single-molecule dynamics.
- Particle detection and tracking were performed using U-Track2 (MATLAB)
- Custom scripts implemented to perform additional analyses on the U-Track2 output, including:
	-	Diffusion coefficients calculation
	-	Trajectory type classification
	-	Spot intensity quantification  
	  These analyses are described in: [JoVe Protocol, Sorzano COS and Martínez-Muñoz L et al, 2019] (https://www.jove.com/video/59314) (DOI:10.3791/59314)
- This approach was previously applied in: [Molecular Cell, Martínez-Muñoz L. et al, 2018] (https://doi.org/10.1016/j.molcel.2018.02.034) 
- Scripts are included here for convenience to allow reviewers to reproduce the results.

## 3. Dual-TIRF analysis for CXCR4 vesicular tracking (folder: Dual_TIRF_CXCR4_tracking)
- Vesicle localization.
- Particle tracking.
- Classification of vesicle dynamics in T cells.
- Quantitative measurements were exported as CSV files for downstream statistical analysis. 
- Graphs and statistical analyses were generated using GraphPad Prism.


*Each analysis folder contains:
- The corresponding scripts (Python, MATLAB, and ImageJ macros).
- Documentation explaining the workflow.
- Information about generated output files.

# Software Environment

- MATLAB R2021a
- Python 3.8.13
- ImageJ2 (version 2.16.0/1.54p)
- GraphPad Prism 10 (used for final figure generation)

# Reproducibility

All scripts required to generate the quantitative data are included.
The scripts generate txt, CSV or Excel tables used for statistical analysis and final figure preparation.
Final figure formatting was performed in GraphPad Prism 10.

# Figure Correspondence

- T cell contact area: Figure 3D-G
- Single-molecule TIRF: Figures 1C-F, 2, 4 and 7; Supplementary Figures 9 and 10.
- Dual-TIRF CXCR4 vesicles: Figures 5 and 6; Supplementary Figures 6 and 7.

