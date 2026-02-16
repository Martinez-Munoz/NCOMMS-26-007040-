# Dual_TIRF_CXCR4_tracking

This folder contains the scripts used to analyze CXCR4 and Rab4/Rab11 cellular structures dynamics.

# # 1. Video preprocessing, particle detection and tracking
- For 2-color datasets, raw split-field videos were first separated into their CXCR4 and Rab4/Rab11 channels, and aligned. 1-color data didn’t require this step.
- Cells were segmented using Cellpose (Stringer et al., 2021).
- Cellular structures (e.g. endosomes, nanoclusters) of varying sized were detected by running scikit-image’s blob detection algorithm (van der Walt et al., 2014), implemented with a Laplacian of Gaussian (LoG) filter.
- Tracking was performed with TrackPy (Allan et al., 2024; Crocker & Grier, 1996) in 3D (x-y, and radius).
- Trajectories were refined by running in parallel ring detection, OpenCV’s Hough Circle Transform, and tracking.

# # 2. Trajectory analysis
- The final CXCR4 trajectories were:
   - Filtered to include a minimum of points
   - Classified as large or small based on a radius threshold
   - In the case of 2-color data, classified as Rab4/Rab11 positive or negative.
- Individual trajectories were characterized by computing:
   - Instantaneous diffusion coefficient (D1-4).
   - Anomalous exponent (a)
   - Total travelled distance

# # 3. Data export
- Quantitative measurements are exported as .csv files for downstream statistical analysis.
- Graphs and statistical analyses were generated using GraphPad Prism 10.


# Scripts Included

1_Cell_segmentation_and_trackings_mac_large_files.py

Reads raw data and returns trajectories for CXCR4 and Rab4/Rab11.

2_Correlations_and_diffusion_plots.py

Reads CXCR4 and Rab4/Rab11 trajectories and returns analysis on their colocalization and dynamics.

3_Merge_csv_files.py

Merges all output csv files for easier use.


1_Cell_segmentation_and_trackings_mac_large_files_1channel.py

Reads raw data and returns trajectories for CXCR4.

2_Correlations_and_diffusion_plots_1channel.py

Reads CXCR4 trajectories and returns analysis on their dynamics.

3_Merge_csv_files_1channel.py

Merges all output csv files for easier use.

