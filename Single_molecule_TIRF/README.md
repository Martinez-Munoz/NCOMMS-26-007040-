#Single_molecule_TIRF

This folder contains the scripts used to analyze CXCR4 and β1 integrin single-molecule dynamics and organization in T cells.

##Single-molecule TIRF Workflow Summary

##1. Particle detection and tracking
- Performed using U-Track2 (MATLAB).
- Tracks were manually verified frame by frame to ensure correct identification.
- Particles incorrectly detected by U-Track2 were that were manually annotated. 

##2. Custom MATLAB analyses on U-Track2 trajectories.
- Trajectories corresponding to incorrectly identified spots were first exlcuded (manually verified).
  - Custom MATLAB scripts were used to extract the X-Y coordinates of each spot along its trajectory.
- From these coordinates, the following parameters were calculated:
- Diffusion coefficients for each trajectory.
   - Long trajectories were classified using the Moment Scaling Spectrum (MSS) analysis.
   - Spot intensity along each trajectory.
- These analyses were first described in: [JoVe Protocol, Sorzano COS and Martínez-Muñoz L et al, 2019] (https://www.jove.com/video/59314) (DOI:10.3791/59314)

##3. Data export
   - Quantitative measurements are exported as TXT files for downstream statistical analysis.
   - Graphs and statistical analyses were generated using GraphPad Prism 10.

##Scripts Included

readTrajectories.m - Reads U-Track2 output and extracts spot trajectories
separateTrajectoriesByLength.m - Splits trajectories based on their length
AnalyzeSpotIntensities.m - Calculates spot intensities along trajectories
gatherDiffusionAndIntensity.m - Combines diffusion coefficients and intensities into a summary table
gatherTrajectoryClassificationAndIntensity.m - Combines trajectory type classifications and spot intensities

DiffusionScripts/  Contains additional MATLAB scripts to calculate diffusion coefficients for all trajectories

-	calculateDiffusion.m
-	FitFunctionGJB10a.m
-	FitMeanMSDconfinedGJB18.m
-	FitMeanMSDdirectedGJB18.m
-	FitMeanMSDlinearGJB18.m
-	FitMeanMSDtalphaGJB18.m
-	LengthAdjustAndSelect12.m
-	MeanMSD_GJB16.m
-	MeanMSDanalysisGJB18.m
-	msddGJB20.m

MomentScalingSpectrum/ Contains MATLAB scripts to perform MSS analysis and classify trajectories

-	AnalyseAllFilesGJB20_MAIN_MSS.m
-	classifyLongTrajectories.m
-	momentScalingSpectrum.m
-	msddGJB20_MSS.m
-	simulateDirectedTrajectories.m



**Output**

- TXT files containing all calculated parameters
- Tables can be used directly for generating figures reported in the manuscript
