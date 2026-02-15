# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 12:21:23 2024

@author: rpons
"""

import numpy as np
import pandas as pd
from tqdm.notebook import tqdm
import matplotlib.pyplot as plt
import trackpy as tp
from PIL import Image, ImageDraw
import cv2
import csv
import glob
from skimage.feature import blob_log
from skimage import feature
import random
from dask_image.imread import imread
import os
from cellpose import models, io
from scipy.ndimage import distance_transform_edt
from skimage.segmentation import watershed
from scipy import ndimage as ndi
from skimage.feature import peak_local_max
from scipy.spatial.distance import cdist
from scipy.ndimage import uniform_filter, gaussian_filter
import logging
import warnings
import seaborn as sns
import re
sns.set_palette("colorblind")
plt.rcParams['figure.dpi'] = 200

# Below, the parameters to check before running the code:
############################################################################################

# In the path below, you should specify where all of your data is stored.
folder_name = os.path.normpath(r"/Volumes/CABIMER 5/Spinning disc/20221111 KIN X4")

# Specify the experiment you would like to analyze/plot.
experiment_to_analyze = "Experiment_1"

# Do you want time plots? (it won't affect the output csv)
time_plots = True # (True=yes, False=no)

# Do you want to plot the alpha's and D's for every video? (e.g. to see variability) (inherited from Lukas's code) (it won't affect the output csv)
alpha_plots = True # (True=yes, False=no)
diffusion_plots = True # (True=yes, False=no)

# Here we will specify the threshold between big spots and small spots.
# You can either select 1 threshold (small vs big), or 2 threhsolds (excludes data between small and big spots)
number_of_radius_thresholds = True # True = 1 threshold, False = 2 thresholds
# Select the threshold in case you set True
radius_threshold = 0.2 #(um)
# Select both threshold in case you set False
lower_radius_threshold = 0.3 #(um)
higher_radius_threshold = 0.6 #(um)

# Specify pixel size and frame rate
pixel_size = 0.0645 # (um)
dt = 0.1 #(s) frame time

#Do you want to export .csv files with all the outputs?
csv_export = True

############################################################################################


## Plotting options
whiskrange = [10, 90] # plot whiskers
rota = -45. # label rotation

condition_to_analyze = None   # Set to "KO" or None for all
time_point_to_analyze = None  # Set to "10 min" or None for all
video_to_analyze = None       # Set to "video1" or None for all

video_folder_paths = []
# Walk through the directory structure
for root, dirs, files in os.walk(folder_name):
    # Extract directory components relative to the root directory
    relative_path = os.path.relpath(root, folder_name)
    path_parts = relative_path.split(os.sep)

    # Apply filters for each layer
    if len(path_parts) >= 1 and experiment_to_analyze and path_parts[0] != experiment_to_analyze:
        continue
    if len(path_parts) >= 2 and condition_to_analyze and path_parts[1] != condition_to_analyze:
        continue
    if len(path_parts) >= 3 and time_point_to_analyze and path_parts[2] != time_point_to_analyze:
        continue

    # Detect "video" level (fourth layer) and stop there
    if len(path_parts) == 4:
        video_folder_name = path_parts[3]
        if video_to_analyze and not video_folder_name.startswith(video_to_analyze):
            continue
        video_folder_paths.append(root)

def import_image(tif_file):
    frames = imread(tif_file)
    video = np.asarray(frames)
    return video
        
def import_tracks(csv_file):
    ## Datatypes used for reading the csv values
    dtypes = {'LABEL': np.str_, 'ID': np.int32, 'TRACK_ID': np.int32, \
            'QUALITY': np.float64, 'POSITION_X': np.float64, \
            'POSITION_Y': np.float64, 'POSITION_Z': np.float64, \
            'POSITION_T': np.float64, 'FRAME': np.int32, \
            'RADIUS': np.float64, 'VISIBILITY': np.int32, \
            'MANUAL_SPOT_COLOR': np.str_, 'MEAN_INTENSITY_CH1': np.float64, \
            'MEDIAN_INTENSITY_CH1': np.float64, 'MIN_INTENSITY_CH1': np.float64, \
            'MAX_INTENSITY_CH1': np.float64, 'TOTAL_INTENSITY_CH1': np.float64, \
            'STD_INTENSITY_CH1': np.float64, 'CONTRAST_CH1': np.float64, \
            'SNR_CH1': np.float64}
        
    #spots = pd.read_csv(os.path.join(os.getcwd(), csv_file), skiprows = [0,1,2], dtype = dtypes, encoding_errors='backslashreplace')
    spots = pd.read_csv(csv_file, skiprows=[1,2,3], dtype=dtypes, encoding_errors='backslashreplace')
    spots.sort_values(by = ['TRACK_ID','POSITION_T'], inplace = True, ignore_index = True)
    return spots

def filter_by_traj_length(trajectories, min_length):
    for idx in trajectories.TRACK_ID.unique():
        if len(trajectories[trajectories["TRACK_ID"] == idx]) < 12:
            trajectories = trajectories[trajectories["TRACK_ID"] != idx]
    return trajectories

def filter_by_colocalization(trajectories_cxcr4, trajectories_rab):
    trajectories_cxcr4["Rab colocalization"] = None
    for idx in trajectories_cxcr4.TRACK_ID.unique():
        track_cxcr4 = trajectories_cxcr4[trajectories_cxcr4["TRACK_ID"] == idx].copy()
        track_cxcr4["RANGE"] = None
        for frame in track_cxcr4.FRAME.unique():
            idxx = int(track_cxcr4.index[track_cxcr4['FRAME'] == frame].to_numpy()[0])
            coord_cxcr4 = track_cxcr4[["POSITION_X", "POSITION_Y"]][track_cxcr4["FRAME"] == frame].to_numpy()
            coord_rab = trajectories_rab[["POSITION_X", "POSITION_Y"]][trajectories_rab["FRAME"] == frame].to_numpy()
            radii_rab = trajectories_rab["RADIUS"][trajectories_rab["FRAME"] == frame].to_numpy()
            distances = np.sqrt(np.sum((coord_rab - coord_cxcr4) ** 2, axis=1))
            ranges = distances-radii_rab
            if (ranges <= 0).any():
                track_cxcr4.at[idxx, "RANGE"] = -1
                trajectories_cxcr4.loc[(trajectories_cxcr4['TRACK_ID'] == idx) & (trajectories_cxcr4['FRAME'] == frame), "Rab colocalization"] = "positive"
            else:
                track_cxcr4.at[idxx, "RANGE"] = 1
                trajectories_cxcr4.loc[(trajectories_cxcr4['TRACK_ID'] == idx) & (trajectories_cxcr4['FRAME'] == frame), "Rab colocalization"] = "negative"
        # if np.sum(track_cxcr4.RANGE) > 1:
        #     trajectories_cxcr4 = trajectories_cxcr4[trajectories_cxcr4["TRACK_ID"] != idx]
        #     #trajectories_cxcr4.at[idx, "Rab colocalization"] = "Positive"
    return trajectories_cxcr4
                
def classify_by_cell(trajectories, cell_mask):
    trajectories["CELL_ID"] = None
    trajectories["CELL_AREA"] = None
    cell_areas = np.empty((0,2))
    for i in np.unique(cell_mask[0]):
        if i != 0:
            cell_area = (cell_mask[0] == i).sum()*pixel_size**2
            cell_areas = np.vstack((cell_areas, [int(str(i)[0]), cell_area]))
    for idx in trajectories.TRACK_ID.unique():
        track_cxcr4 = trajectories[trajectories["TRACK_ID"] == idx].copy()
        x_avg = int(np.mean(track_cxcr4.POSITION_X)/pixel_size)
        y_avg = int(np.mean(track_cxcr4.POSITION_Y)/pixel_size)
        if int(str(cell_mask[0, y_avg, x_avg])[0]) != 0:
            trajectories.loc[trajectories['TRACK_ID'] == idx, 'CELL_ID'] = int(str(cell_mask[0, y_avg, x_avg])[0])
            trajectories.loc[trajectories['TRACK_ID'] == idx, 'CELL_AREA'] = cell_areas[:,1][cell_areas[:,0]==int(str(cell_mask[0, y_avg, x_avg])[0])][0]
        else:
            trajectories = trajectories[trajectories["TRACK_ID"] != idx]
    return trajectories

def manage_missed_localizations(df, interpolation):
    filled_data = []
    
    for track_id, group in df.groupby('TRACK_ID'):
        group = group.sort_values('FRAME')  # Ensure sorted by FRAME
        frames = group['FRAME'].tolist()
        full_range = range(frames[0], frames[-1] + 1)  # Complete frame range
        
        # Add existing rows
        filled_data.append(group)
        
        # Find and fill missing frames
        for frame in full_range:
            if frame not in frames:
                # Get the rows before and after the missing frame
                before = group[group['FRAME'] < frame].iloc[-1]
                after = group[group['FRAME'] > frame].iloc[0]
                
                if interpolation:
                    # Compute midpoint
                    midpoint = {
                        'TRACK_ID': track_id,
                        'FRAME': frame,
                        'POSITION_X': (before['POSITION_X'] + after['POSITION_X']) / 2,
                        'POSITION_Y': (before['POSITION_Y'] + after['POSITION_Y']) / 2
                    }
                else:
                    midpoint = {
                        'TRACK_ID': track_id,
                        'FRAME': frame,
                        'POSITION_X': np.nan,
                        'POSITION_Y': np.nan
                    }
                filled_data.append(pd.DataFrame([midpoint]))
    
    # Combine all data and sort
    filled_df = pd.concat(filled_data).sort_values(['TRACK_ID', 'FRAME']).reset_index(drop=True)
    return filled_df

def getsampletype(s):
    if ('CTL' in s):
        return 'CTL'
    if ('SDF' in s):
        return 'CTL'
    if ('KO' in s):
        return 'KO'
    if ('DOX' in s):
        return 'DOX'
    else:
        return s.split('_')[0]

def calc_D_and_alpha (traj, px_size = 1, dt = 1, plotmsd = False):
    
    if traj.empty:
        return
    import trackpy as tpy

    max_lt_alpha = int(round(len(traj)*0.3)) # lagtime used for alpha calc. , 30 % of trace
    traj_msd = tpy.motion.msd(traj, px_size, 1./dt, max_lagtime = max_lt_alpha)['msd']
    traj_msd.index = dt*traj_msd.index

    
    if plotmsd:
        plt.figure()
        plt.ylabel(r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]')
        plt.xlabel('lag time $t$');

    fit_parameters = tpy.utils.fit_powerlaw(traj_msd, plot = plotmsd) # from documentation of trackpy, MSD = Dt**alpha
    
    
    alpha = fit_parameters['n'].iloc[0]
    D = fit_parameters['A'].iloc[0]/4.
    
    index = ['D_longterm','alpha']
    
    return pd.Series((D, alpha), index = index)

def calc_D_1_4 (traj, px_size = 1, dt = 1, plotmsd = False):
    
    if traj.empty:
        return
    import trackpy as tpy
    
    max_lt_alpha = 4 # lagtime used for alpha calc. , 1-4 for D_1_4
    traj_msd = tpy.motion.msd(traj, px_size, 1./dt, max_lagtime = max_lt_alpha)['msd']
    traj_msd.index = dt*traj_msd.index
    
    fit_coeff = np.polynomial.polynomial.polyfit(x = traj_msd.index, y = traj_msd, deg = 1) ## Linear fit to alpha
    
    x_0 = fit_coeff[0] ## y intercept
    D_1_4 = fit_coeff[1]/4. ## Slope of the MSD 
    
    
    
    if D_1_4 < 10**(-5):
        #print("D_1_4: Stationary trajectory, D < 10**(-5.)")
        D_1_4 = 10**(-5)
    
    if plotmsd:
        
        y_fit = (4. * D_1_4 * traj_msd.index) + x_0
        
        from matplotlib import pyplot as plt
        plt.figure()
        plt.ylabel(r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]')
        plt.xlabel('lag time $t$');
        
        
        plt.plot(traj_msd.index, traj_msd, 'ro', traj_msd.index, y_fit, 'b-')
        plt.show()
        
    index = ['D_1_4','D_1_4_intercept']
    return pd.Series((D_1_4, x_0), index = index)

def calc_D_1_4_fix_intercept (traj, px_size = 1, dt = 1, plotmsd = False):
    
    # The intercept is fixed to the range of the medians observed in all samples. Big -> very low. Small -> Around 0.05
    
    if traj.empty:
        return
    import trackpy as tpy
    
    max_lt_alpha = 4 # lagtime used for alpha calc. , 1-4 for D_1_4
    traj_msd = tpy.motion.msd(traj, px_size, 1./dt, max_lagtime = max_lt_alpha)['msd']
    traj_msd.index = dt*traj_msd.index
    
    def func(x, a, b):
        return a + b * x
    
    from scipy.optimize import curve_fit
    fit_coeff, _ = curve_fit(func, traj_msd.index, traj_msd, bounds = ([0., 10**(-10)], [0.05, 10**(10)]))  # bounds=([a_lower,b_lower],[a_upper,b_upper])

    x_0 = fit_coeff[0] ## y intercept
    D_1_4 = fit_coeff[1]/4. ## Slope of the MSD 
    
    if plotmsd:
        
        y_fit = (4. * D_1_4 * traj_msd.index) + x_0
        
        from matplotlib import pyplot as plt
        plt.figure()
        plt.ylabel(r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]')
        plt.xlabel('lag time $t$');
        
        
        plt.plot(traj_msd.index, traj_msd, 'ro', traj_msd.index, y_fit, 'b-')
        plt.show()
    
    index = ['D_1_4_fix_intercept']
    return pd.Series((D_1_4), index = index)

def calc_mean_inst_vel (traj, px_size = 1, dt = 1, labels = ['FRAME','POSITION_X','POSITION_Y']):
    
    if traj.empty:
        return
    t, x, y = traj[labels[0]].to_numpy(), traj[labels[1]].to_numpy(), traj[labels[2]].to_numpy()
    
    dx = np.diff(x)
    dy = np.diff(y)
    dt2 = np.diff(t)*dt
    
    dist = np.sqrt(dx**2 + dy**2)
    inst_vel = dist/dt2
    
    mean_inst_vel = np.mean(inst_vel)
    dist_travelled = np.sum(dist)
    se_dist = np.sqrt((x[-1]-x[0])**2 + (y[-1]-y[0])**2)
    
    index = ['mean_inst_vel','travelled_dist', 'start_end_dist']
    
    return pd.Series((mean_inst_vel, dist_travelled, se_dist), index = index)


def powlaw(x, a, b) :
    return a * np.power(x, b)
def linlaw(x, a, b) :
    return a + x * b

def curve_fit_log(xdata, ydata, p0 = [0.1, 0.5]) :
    from scipy.optimize import curve_fit
    """Fit data to a power law with weights according to a log scale"""
    # Weights according to a log scale
    # Apply fscalex
    xdata_log = np.log10(xdata)
    # Apply fscaley
    ydata_log = np.log10(ydata)
    # Fit linear
    popt_log, pcov_log = curve_fit(linlaw, xdata_log, ydata_log, p0 = p0)
    #print(popt_log, pcov_log)
    # Apply fscaley^-1 to fitted data
    ydatafit_log = np.power(10, linlaw(xdata_log, *popt_log))
    # There is no need to apply fscalex^-1 as original data is already available
    return (popt_log, pcov_log, ydatafit_log)

        
def calc_MSS_slope (traj, labels = ['FRAME','POSITION_X','POSITION_Y']):
    
    if traj.empty:
        return
    
    import numpy as np
    ## CHeck https://stackoverflow.com/questions/41109122/fitting-a-curve-to-a-power-law-distribution-with-curve-fit-does-not-work
    # tdiff = max(traj[labels[0]])-min(traj[labels[0]]) # fill missing values in trace with nan to fill gaps
    # if len(traj[labels[0]]) < tdiff:
    #     traj.loc[traj.shape[0]] = [np.nan]*traj.shape[1]
        
    # traj = traj.sort(by = labels[0]) # sort by time
    
    #import pdb; pdb.set_trace()            
    t, x, y = traj[labels[0]].to_numpy(), traj[labels[1]].to_numpy(), traj[labels[2]].to_numpy()
    
    
    
    # Vega et al. 2018 Biopys J. Multistep Track Segmentation and Motion Classification for Transient Mobility Analysis
    max_lt = min(np.floor(0.5 + 0.25*len(traj.index)).astype(int), 30) # lagtime used for MSS calc., 25 % of trace or max 30 frames
    #max_lt = int(round(len(traj)*0.3)) # lagtime used for alpha calc. , 30 % of trace
    alphas = [0]
    
    for j in np.arange(start = 1, stop = 7): ## Calculate moments 1-6
        #traj_msd = tpy.motion.msd(traj, px_size, 1./dt, max_lagtime = max_lt_alpha)['msd']
        moments_dt = []
        for dt in np.arange(start = 1, stop = max_lt+1):
            moment = 0.
            dts = 0
            for m in np.arange(start = 0, stop = (len(t)-dt)): ## len(t)-1-dt, but np arange does not include the stop value
                #if (t[m+dt]-t[m] == dt): ## There may be gaps in the trace, this checks
                moment += ((x[m+dt]-x[m])**2 + (y[m+dt]-y[m])**2)**(j/2.)
                dts += 1
                    
                    ## Exponentiation by squaring might be more efficient here, but is more tricky with the gaps.
                    # Maybe this can be better recoded with a numpy array [t x y]
                    
            #moment /= float(len(t)-dt)
            moment /= dts
            moments_dt.append(moment)
        #traj_moment = sp.stats.moment(traj, moment = i, axis = 0, nan_policy = 'propagate')
        #scipy moment(n, dfn, dfd, nc, loc = 0, scale = 1)
        
        
        # traj_moments = pd.Series(moments_dt, index = [frametime*x for x in np.arange(start = 1,stop = max_lt+1)])
        # alpha_fit_parameters = tpy.utils.fit_powerlaw(traj_moments, plot = plotmoment) # from documentation of trackpy, MSD = Dt**alpha
        # alpha = alpha_fit_parameters['n'].iloc[0]
        
        popt_log, pcov_log, ydatafit_log = curve_fit_log([x for x in np.arange(start = 1,stop = max_lt+1)], moments_dt)
        alpha = popt_log[1]
        
        alphas.append(alpha)

        # plotmsd = False
        # if plotmsd:
        #     plt.figure()
        #     plt.ylabel(r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]')
        #     plt.xlabel('lag time $t$');
    x = np.arange(start = 0, stop = len(alphas))
    coeff = np.polyfit(x, alphas, 1) ## Linear fit to alpha
    
    return coeff[0]

def calculate_diffusion_parameters(trajectories):
    
    keys = [key for key, _ in trajectories.groupby('TRACK_ID')] 
    grouped = trajectories.groupby('TRACK_ID')
    cells_by_key = {
    key: (group['CELL_ID'].iloc[0], group['CELL_AREA'].iloc[0]) for key, group in grouped
    }

    tracks = pd.DataFrame({
    'TRACK_ID': keys,
    'CELL_ID': [cells_by_key[key][0] for key in keys],
    'CELL_AREA': [cells_by_key[key][1] for key in keys]
    })
    
    tracks['filename'] = relative_name + '_' + tracks['CELL_ID'].astype(str)
    tracks['sampletype'] = relative_path.split(os.sep)[1]
    tracks['time'] = time[0]
    tracks['track_length'] = trajectories.groupby('TRACK_ID').size().reset_index(drop = True)
    
    tpy_spots = trajectories.rename(columns = {'FRAME': 'frame', 'POSITION_X': 'x','POSITION_Y': 'y', 'TRACK_ID':'particle'})
    
    D_alpha = tpy_spots.groupby('particle').apply(
        lambda x: calc_D_and_alpha(x, 1, dt)
        ).reset_index(drop = True)
    
    D_1_4 = tpy_spots.groupby('particle').apply(
        lambda x: calc_D_1_4(x, 1, dt)
        ).reset_index(drop = True)
    
    D_1_4_fi = tpy_spots.groupby('particle').apply(
        lambda x: calc_D_1_4_fix_intercept(x, 1, dt)
        ).reset_index(drop = True)

    vel_dist = tpy_spots.groupby('particle').apply(
        lambda x: calc_mean_inst_vel(x, 1, dt, ['frame', 'x', 'y'])
        ).reset_index(drop = True)
    
    tracks['MSS_slope'] = tpy_spots.groupby('particle').apply(
        lambda x: calc_MSS_slope(x, ['frame', 'x', 'y'])
        ).reset_index(drop = True)
    
    tracks = pd.concat([tracks, D_alpha, D_1_4, D_1_4_fi, vel_dist], axis = 1)
    
    tracks['track_intensity_max'] = trajectories.groupby('TRACK_ID')['MAX_INTENSITY_CH1'].max().reset_index(drop = True)
    tracks['track_intensity_mean'] = trajectories.groupby('TRACK_ID')['MEAN_INTENSITY_CH1'].mean().reset_index(drop = True)
    tracks['radius_mean'] = trajectories.groupby('TRACK_ID')['RADIUS'].mean().reset_index(drop = True)

    return tracks

# def sorting_key(label):
#     # Check if the label contains numbers
#     if re.search(r"\d", label):
#         # Extract numeric part and return as sorting key (ensure it's an integer)
#         numeric_part = int("".join(filter(str.isdigit, label)))
#         return (2, numeric_part)  # Group numeric labels after non-numeric ones
#     else:
#         # Check for specific keywords in non-numeric labels
#         if "CTL" in label:
#             return (0, 0)  # CTL comes first
#         elif "KO" in label:
#             return (0, 1)  # KO comes next
#         else:
#             return (0, 2)  # Any other non-numeric labels come last

def sorting_key(label):
    parts = label.split()  # Split into parts
    group = 0 if parts[0] == 'CTL' else 1  # CTL comes before KO
    size = 0 if parts[1] == 'small' else 1  # small comes before big
    return (group, size)


all_outputs = []
print(f"Plotting...")
for i in range(len(video_folder_paths)):
    # import pdb
    # pdb.set_trace()
    relative_path = os.path.relpath(video_folder_paths[i], folder_name)
    relative_experiment = relative_path.split(os.sep)[0]
    relative_name = relative_path.split(os.sep)[0] + "_" + relative_path.split(os.sep)[1] + "_" + relative_path.split(os.sep)[2] + "_" + relative_path.split(os.sep)[3]
    save_path = os.path.join(folder_name, relative_experiment, relative_name)
    print(f"({i+1}/{len(video_folder_paths)}) Analyzing video {relative_path}...")
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning, module="pims.api")
        try:
            if os.path.exists(save_path + '_Cellpose_cell_mask.tif') and os.path.exists(save_path + '_Tracks_CXCR4.csv'):
                time = relative_path.split(os.sep)[2]
                time = int(re.search(r'\d+', time).group())
                time = np.array([time])
                
                #Importing the cell mask and tracks
                cell_mask = import_image(save_path + '_Cellpose_cell_mask.tif')
                trajectories_cxcr4 = import_tracks(save_path + '_Tracks_CXCR4.csv')
                #trajectories_rab = import_tracks(save_path + '_Tracks_Rab.csv')
            
                #Filtering CXCR4 tracks by trajectory length and by colocalization with Rab trajectories
                trajectories_cxcr4 = filter_by_traj_length(trajectories_cxcr4, min_length=12)
                #trajectories_rab = filter_by_traj_length(trajectories_rab, min_length=12)
                
                #trajectories_cxcr4 = filter_by_colocalization(trajectories_cxcr4, trajectories_rab)
                
                #Classifying trajectories based on the cell mask (= it tells to which cell the tracks belong to, in case there are multiple cells in the same video)
                trajectories_cxcr4 = classify_by_cell(trajectories_cxcr4, cell_mask)
                
                if len(trajectories_cxcr4) != 0:
                    #Calculating diffusion parameters (D1-4, alpha, D longterm, velocity...)
                    output = calculate_diffusion_parameters(trajectories_cxcr4)
                    if number_of_radius_thresholds:
                        output['sampletype'] = output['sampletype'] + output['radius_mean'].apply(lambda size: ' big' if size > radius_threshold else ' small')
                    if not number_of_radius_thresholds:
                        output['sampletype'] = output['sampletype'] + output['radius_mean'].apply(lambda size: ' big' if size > higher_radius_threshold else ' small')
                        output = output[(output["radius_mean"] <= lower_radius_threshold) | (output["radius_mean"] >= higher_radius_threshold)]
                    output = output.sort_values('sampletype')
                    if csv_export:
                        relative_experiment = relative_path.split(os.sep)[0]
                        relative_name = relative_path.split(os.sep)[0] + "_" + relative_path.split(os.sep)[1] + "_" + relative_path.split(os.sep)[2] + "_" + relative_path.split(os.sep)[3]
                        save_path = os.path.join(folder_name, relative_experiment, relative_name)
                        for cell_id, group in output.groupby('CELL_ID'):
                            group.to_csv(save_path + f'_diffusion_parameters_cell_{cell_id}.csv', encoding='utf-8', index=False)
                        
                    all_outputs.append(output)
                else:
                    print("No CXCR4 trajectories for this cell. Possibly because no cell was segmented.")
        except Exception as e:
            # Handle the error and skip the iteration
            print(f"Error processing data from {relative_name}. Skipping to the next video.")
all_outputs = pd.concat(all_outputs, ignore_index=True)
#all_outputs.to_csv(os.path.join(folder_name, relative_experiment, relative_experiment) + f'_diffusion_parameters_combined_x.csv', encoding='utf-8', index=False)


#all_outputs["sort_key"] = all_outputs["sampletype"].apply(sorting_key)  # Apply sorting key function
#all_outputs = all_outputs.sort_values("sort_key").drop(columns="sort_key").reset_index(drop=True)  # Sort and cleanup

all_outputs_s = all_outputs[all_outputs["radius_mean"] < radius_threshold].copy()
all_outputs_b = all_outputs[all_outputs["radius_mean"] >= radius_threshold].copy()

all_outputs = all_outputs.sort_values(by='sampletype', key=lambda col: col.map(sorting_key))
all_outputs_s = all_outputs_s.sort_values(by='sampletype', key=lambda col: col.map(sorting_key))
all_outputs_b = all_outputs_b.sort_values(by='sampletype', key=lambda col: col.map(sorting_key))




if alpha_plots:
    
    fig, ax = plt.subplots()
    sns.boxplot(x = 'sampletype', y = 'alpha', data = all_outputs,
                whis = whiskrange, width = .6, hue='sampletype', palette = 'vlag', showfliers = False)
    p = sns.stripplot(x = 'sampletype', y = 'alpha', data = all_outputs, hue='sampletype', palette = 'vlag',
                  size = 4, color = '.3', linewidth = 0.5)
    p.axhline(1)
    #plt.xticks(rotation = rota)
    plt.subplots_adjust(bottom = 0.3)
    plt.setp( ax.xaxis.get_majorticklabels(), rotation = rota, ha = 'left', rotation_mode = 'anchor')
    p.set(xlabel = 'Sample', ylabel = r'$\alpha$')
    sns.despine()
    

    fig, ax = plt.subplots()
    sns.boxplot(x = 'sampletype', y = 'D_1_4_fix_intercept', data = all_outputs,
                whis = whiskrange, width = .6, hue='sampletype', palette = 'vlag', showfliers = False)
    p = sns.stripplot(x = 'sampletype', y = 'D_1_4_fix_intercept', data = all_outputs, hue='sampletype', palette = 'vlag',
                  size = 4, color = '.3', linewidth = 0.5, )
    p.set_yscale('log')
    #plt.xticks(rotation = rota)
    plt.subplots_adjust(bottom = 0.3)
    plt.setp( ax.xaxis.get_majorticklabels(), rotation = rota, ha = 'left', rotation_mode = 'anchor')
    p.set(xlabel = 'Sample', ylabel = r'D_1_4_fix_intercept')
    sns.despine()
    plt.close(fig)
        
if diffusion_plots: 
       
    ### Plot D per file to see variability
    fig, ax = plt.subplots()
    sns.boxplot(x = 'filename', y = 'D_1_4', data = all_outputs_b,
                whis = whiskrange, width = .6, hue='sampletype', palette = 'vlag', showfliers = False)
    p = sns.stripplot(x = 'filename', y = 'D_1_4', data = all_outputs_b, hue='sampletype', palette = 'vlag',
                  size = 4, color = '.3', linewidth = 0.5)
    #plt.xticks(rotation = rota)
    plt.subplots_adjust(left = 0.2, bottom = 0.3)
    p.set_yscale('log')
    plt.setp( ax.xaxis.get_majorticklabels(), rotation = rota, ha = 'left', rotation_mode = 'anchor')
    p.set(xlabel = 'Sample', ylabel = r'$D_{1-4}\ (um^2/s)$')
    sns.despine()
    
    fig, ax = plt.subplots()
    sns.boxplot(x = 'filename', y = 'D_1_4_fix_intercept', data = all_outputs_s,
                whis = whiskrange, width = .6, hue='sampletype', palette = 'vlag', showfliers = False)
    p = sns.stripplot(x = 'filename', y = 'D_1_4', data = all_outputs_s, hue='sampletype', palette = 'vlag',
                  size = 4, color = '.3', linewidth = 0.5)
    #plt.xticks(rotation = rota)
    plt.subplots_adjust(left = 0.2, bottom = 0.3)
    p.set_yscale('log')
    plt.setp( ax.xaxis.get_majorticklabels(), rotation = rota, ha = 'left', rotation_mode = 'anchor')
    p.set(xlabel = 'Sample', ylabel = r'$D_{1-4}\ (um^2/s)$')
    sns.despine()
    
    fig, ax = plt.subplots()
    sns.boxplot(x = 'filename', y = 'alpha', data = all_outputs_b,
                whis = whiskrange, width = .6, hue='sampletype', palette = 'vlag', showfliers = False)
    p = sns.stripplot(x = 'filename', y = 'alpha', data = all_outputs_b, hue='sampletype', palette = 'vlag',
                  size = 4, color = '.3', linewidth = 0.5)
    p.axhline(1)
    #plt.xticks(rotation = rota)
    plt.subplots_adjust(bottom = 0.3)
    plt.setp( ax.xaxis.get_majorticklabels(), rotation = rota, ha = 'left', rotation_mode = 'anchor')
    p.set(xlabel = 'Sample', ylabel = r'$\alpha$')
    sns.despine()
    
    fig, ax = plt.subplots()
    sns.boxplot(x = 'filename', y = 'alpha', data = all_outputs_s,
                whis = whiskrange, width = .6, hue='sampletype', palette = 'vlag', showfliers = False)
    p = sns.stripplot(x = 'filename', y = 'alpha', data = all_outputs_s, hue='sampletype', palette = 'vlag',
                  size = 4, color = '.3', linewidth = 0.5)
    p.axhline(1)
    #plt.xticks(rotation = rota)
    plt.subplots_adjust(bottom = 0.3)
    plt.setp( ax.xaxis.get_majorticklabels(), rotation = rota, ha = 'left', rotation_mode = 'anchor')
    p.set(xlabel = 'Sample', ylabel = r'$\alpha$')
    sns.despine()
    
    
    ### Plot D, alpha and MSS slope per sampletype
    
    fig, ax = plt.subplots()
    sns.boxplot(x = 'sampletype', y = 'D_1_4', data = all_outputs,
                whis = whiskrange, width = .6, hue='sampletype', palette = 'vlag', showfliers = False)
    p = sns.stripplot(x = 'sampletype', y = 'D_1_4', data = all_outputs, hue='sampletype', palette = 'vlag',
                  size = 4, color = '.3', linewidth = 0.5)
    #plt.xticks(rotation = rota)
    plt.subplots_adjust(left = 0.2, bottom = 0.3)
    p.set_yscale('log')
    plt.setp( ax.xaxis.get_majorticklabels(), rotation = rota, ha = 'left', rotation_mode = 'anchor')
    p.set(xlabel = 'Sample', ylabel = r'$D_{1-4}\ (um^2/s)$')
    sns.despine()
    
    fig, ax = plt.subplots()
    sns.boxplot(x = 'sampletype', y = 'D_longterm', data = all_outputs,
                whis = whiskrange, width = .6, hue='sampletype', palette = 'vlag', showfliers = False)
    p = sns.stripplot(x = 'sampletype', y = 'D_longterm', data = all_outputs, hue='sampletype', palette = 'vlag',
                  size = 4, color = '.3', linewidth = 0.5)
    #plt.xticks(rotation = rota)
    plt.subplots_adjust(left = 0.2, bottom = 0.3)
    p.set_yscale('log')
    plt.setp( ax.xaxis.get_majorticklabels(), rotation = rota, ha = 'left', rotation_mode = 'anchor')
    p.set(xlabel = 'Sample', ylabel = r'D longterm $(um^2/s^{\alpha})$')
    sns.despine()
    
    fig, ax = plt.subplots()
    sns.boxplot(x = 'sampletype', y = 'alpha', data = all_outputs,
                whis = whiskrange, width = .6, hue='sampletype', palette = 'vlag', showfliers = False)
    p = sns.stripplot(x = 'sampletype', y = 'alpha', data = all_outputs, hue='sampletype', palette = 'vlag',
                  size = 4, color = '.3', linewidth = 0.5, )
    p.axhline(1)
    #plt.xticks(rotation = rota)
    plt.subplots_adjust(bottom = 0.3)
    plt.setp( ax.xaxis.get_majorticklabels(), rotation = rota, ha = 'left', rotation_mode = 'anchor')
    p.set(xlabel = 'Sample', ylabel = r'$\alpha$')
    sns.despine()
    
    fig, ax = plt.subplots()
    sns.boxplot(x = 'sampletype', y = 'MSS_slope', data = all_outputs,
                whis = whiskrange, width = .6, hue='sampletype', palette = 'vlag', showfliers = False)
    p = sns.stripplot(x = 'sampletype', y = 'MSS_slope', data = all_outputs, hue='sampletype', palette = 'vlag',
                  size = 4, color = '.3', linewidth = 0.5, )
    p.axhline(0.5)
    #plt.xticks(rotation = rota)
    plt.subplots_adjust(bottom = 0.3)
    plt.setp( ax.xaxis.get_majorticklabels(), rotation = rota, ha = 'left', rotation_mode = 'anchor')
    p.set(xlabel = 'Sample', ylabel = r'MSS slope')
    sns.despine()


import numpy as np
bins = np.logspace(-4, 0, num = 21)
fig, ax = plt.subplots()
p = sns.histplot(x = 'D_1_4', data = all_outputs[all_outputs['sampletype'].isin(['KO small', 'CTL small'])], hue = 'sampletype', element = 'step',
                 stat = 'percent', common_norm = False, \
                     multiple = 'layer', shrink = 1, bins = bins)
ax.set_xscale('log')
p.set(xlabel = r'$D_{1-4}\ (um^2/s)$')

fig, ax = plt.subplots()
p = sns.histplot(x = 'D_1_4', data = all_outputs[all_outputs['sampletype'].isin(['KO big', 'CTL big'])], hue = 'sampletype', element = 'step',
                 stat = 'percent', common_norm = False, \
                     multiple = 'layer', shrink = 1, bins = bins)
ax.set_xscale('log')
p.set(xlabel = r'$D_{1-4}\ (um^2/s)$')
    
bins = np.linspace(-0.2, 2.0, num = 23)
fig, ax = plt.subplots()
p = sns.histplot(x = 'alpha', data = all_outputs[all_outputs['sampletype'].isin(['KO small', 'CTL small'])], hue = 'sampletype', element = 'step',
                 stat = 'percent', common_norm = False, \
                     multiple = 'layer', shrink = 1, bins = bins)
p.set(xlabel = r'$\alpha$')
fig, ax = plt.subplots()
p = sns.histplot(x = 'alpha', data = all_outputs[all_outputs['sampletype'].isin(['KO big', 'CTL big'])], hue = 'sampletype', element = 'step',
                 stat = 'percent', common_norm = False, \
                     multiple = 'layer', shrink = 1, bins = bins)
p.set(xlabel = r'$\alpha$')


    
### Histograms for the diffusion coefficients 1_4 and alpha

import numpy as np
bins = np.logspace(-4, 0, num = 21)
fig, ax = plt.subplots()
p = sns.histplot(x = 'D_longterm', data = all_outputs[all_outputs['sampletype'].isin(['KO big', 'CTL big'])], hue = 'sampletype', element = 'step',
                 stat = 'percent', common_norm = False, \
                     multiple = 'layer', shrink = 1, bins = bins)
ax.set_xscale('log')
p.set(xlabel = r'D longterm $(um^2/s^{\alpha})$')

fig, ax = plt.subplots()
p = sns.histplot(x = 'D_longterm', data = all_outputs[all_outputs['sampletype'].isin(['KO small', 'CTL small'])], hue = 'sampletype', element = 'step',
                 stat = 'percent', common_norm = False, \
                     multiple = 'layer', shrink = 1, bins = bins)
ax.set_xscale('log')
p.set(xlabel = r'D longterm $(um^2/s^{\alpha})$')
    


# Plot alpha vs kdeplot+scatterplot 

y = 'D_longterm'
x = 'alpha'

fig, ax = plt.subplots(1,2, sharex = 'all', sharey = 'all')
p = sns.scatterplot(y = y, x = x, data = all_outputs_s, ax = ax[0],
              hue = 'sampletype')
p = sns.kdeplot(y = y, x = x, data = all_outputs_s, ax = ax[0],
              hue = 'sampletype', log_scale = (False, True), \
                  levels = 3, common_norm = False)
p.set(xlabel = r'$\alpha$', ylabel = r'D longterm $(um^2/s^{\alpha})$')
p = sns.scatterplot(y = y, x = x, data = all_outputs_b, ax = ax[1],
              hue = 'sampletype')
p = sns.kdeplot(y = y, x = x, data = all_outputs_b, ax = ax[1],
              hue = 'sampletype', log_scale = (False, True), \
                  levels = 3, common_norm = False)
plt.tight_layout()
p.set(xlabel = r'$\alpha$', ylabel = r'D longterm $(um^2/s^{\alpha})$')
sns.despine()



# Plot alpha and D vs brightness kdeplot+scatterplot 

fig, ax = plt.subplots(1,2, sharex = 'all') # sharey = 'all')
p = sns.kdeplot(y = 'track_intensity_max', x = 'alpha', data = all_outputs_s, ax = ax[0],
              hue = 'sampletype', log_scale = (False, False), \
                  levels = 3, common_norm = False)
p = sns.scatterplot(y = 'track_intensity_max', x = 'alpha', data = all_outputs_s, ax = ax[0],
              hue = 'sampletype')
p.set(xlabel = r'$\alpha$', ylabel = r'Max track intensity (arb. units)')
p = sns.kdeplot(y = 'track_intensity_max', x = 'alpha', data = all_outputs_b, ax = ax[1],
              hue = 'sampletype', log_scale = (False, False), \
                  levels = 3, common_norm = False)
p = sns.scatterplot(y = 'track_intensity_max', x = 'alpha', data = all_outputs_b, ax = ax[1],
              hue = 'sampletype')
plt.tight_layout()
p.set(xlabel = r'$\alpha$', ylabel = r'Max track intensity (arb. units)')
sns.despine()

fig, ax = plt.subplots(1,1, sharex = 'all', sharey = 'all')
p = sns.scatterplot(y = 'track_intensity_max', x = 'D_longterm', data = all_outputs, ax = ax,
              hue = 'sampletype')
p = sns.kdeplot(y = 'track_intensity_max', x = 'D_longterm', data = all_outputs, ax = ax,
              hue = 'sampletype', log_scale = (True, False), \
                  levels = 3, common_norm = False)
p.set(xlabel = r'D longterm $(um^2/s^{\alpha})$', ylabel = r'Max track intensity (arb. units)')
plt.tight_layout()
sns.despine()


#Plot alpha and D vs mean track radius (scatters)
data = all_outputs.copy()
data['sampletype'] = data['sampletype'].str.replace(r' big| small', '', regex=True)

fig, ax = plt.subplots(1,1, sharex = 'all') # sharey = 'all')

p = sns.kdeplot(y = 'radius_mean', x = 'alpha', data = data, ax = ax, hue = 'sampletype',
              log_scale = (False, False), \
                  levels = 3, common_norm = False)
p = sns.scatterplot(y = 'radius_mean', x = 'alpha', data = data, ax = ax, hue = 'sampletype', alpha=0.5
              )
p.set(xlabel = r'$\alpha$', ylabel = r'Mean track radius (um)')
p.set(xlabel = r'$\alpha$', ylabel = r'Mean track radius (um)')
sns.despine()

fig, ax = plt.subplots(1,1, sharex = 'all', sharey = 'all')
p = sns.scatterplot(y = 'radius_mean', x = 'D_longterm', data = data, ax = ax, hue = 'sampletype', alpha=0.5
              )
p = sns.kdeplot(y = 'radius_mean', x = 'D_longterm', data = data, ax = ax, hue = 'sampletype',
              log_scale = (True, False), \
                  levels = 3, common_norm = False)
p.set(xlabel = r'D longterm $(um^2/s^{\alpha})$', ylabel = r'Mean track radius (um)')
plt.tight_layout()
sns.despine()



#Plot size distribution
fig, ax = plt.subplots(1,1, sharex = 'all') # sharey = 'all')
p = sns.histplot(x = 'radius_mean', data = data, hue = 'sampletype', element = 'step',
                 stat = 'percent', common_norm = False, \
                     multiple = 'layer', shrink = 1, bins = 10)
p.set(xlabel = r'Mean particle radius per track (um)', ylabel = r'Percent')
sns.despine()



#Plot size densities
grouped_s = all_outputs_s.groupby("filename").size().reset_index(name="number_spots")
grouped_s["sampletype"] = grouped_s["filename"].map(all_outputs_s.groupby("filename")["sampletype"].first())
grouped_s["CELL_AREA"] = grouped_s["filename"].map(all_outputs_s.groupby("filename")["CELL_AREA"].first())
grouped_s["density"] = grouped_s["number_spots"] / grouped_s["CELL_AREA"]

grouped_b = all_outputs_b.groupby("filename").size().reset_index(name="number_spots")
grouped_b["sampletype"] = grouped_b["filename"].map(all_outputs_b.groupby("filename")["sampletype"].first())
grouped_b["CELL_AREA"] = grouped_b["filename"].map(all_outputs_b.groupby("filename")["CELL_AREA"].first())
grouped_b["density"] = grouped_b["number_spots"] / grouped_b["CELL_AREA"]

grouped = pd.concat([grouped_s, grouped_b], ignore_index=True)

fig, ax = plt.subplots()
sns.boxplot(
    x="sampletype", y="density", data=grouped, 
    whis=whiskrange, width=0.6, hue="sampletype", 
    palette="vlag", showfliers=False, ax=ax)
sns.stripplot(
    x="sampletype", y="density", data=grouped, 
    hue="sampletype", palette="vlag", size=4, 
    color=".3", linewidth=0.5, dodge=True, ax=ax)
plt.subplots_adjust(bottom=0.3)
plt.setp(ax.xaxis.get_majorticklabels(), rotation=rota, ha='left', rotation_mode='anchor')
ax.set(xlabel="Sample", ylabel=r'Density $(um^{-2})$')
sns.despine()

    
if time_plots:
    all_outputs["sampletype"] = all_outputs.apply(
    lambda row: row['sampletype'] + " " + str(row['time']) + "min", axis=1)
    all_outputs_s["sampletype"] = all_outputs_s.apply(
    lambda row: row['sampletype'] + " " + str(row['time']) + "min", axis=1)
    all_outputs_b["sampletype"] = all_outputs_b.apply(
    lambda row: row['sampletype'] + " " + str(row['time']) + "min", axis=1)
    
    all_outputs = all_outputs.sort_values(by='time')
    all_outputs_s = all_outputs_s.sort_values(by='time')
    all_outputs_b = all_outputs_b.sort_values(by='time')
    
    bins = np.logspace(-4, 0, num = 21)
    fig, ax = plt.subplots()
    data = all_outputs[all_outputs['sampletype'].str.contains('KO small', case=False, na=False)]
    p = sns.histplot(x = 'D_1_4', data = data, hue = 'sampletype', element = 'step',
                     stat = 'percent', common_norm = False, \
                         multiple = 'layer', shrink = 1, bins = bins)
    ax.set_xscale('log')
    p.set(xlabel = r'$D_{1-4}\ (um^2/s)$')

    fig, ax = plt.subplots()
    data = all_outputs[all_outputs['sampletype'].str.contains('KO big', case=False, na=False)]
    p = sns.histplot(x = 'D_1_4', data = data, hue = 'sampletype', element = 'step',
                     stat = 'percent', common_norm = False, \
                         multiple = 'layer', shrink = 1, bins = bins)
    ax.set_xscale('log')
    p.set(xlabel = r'$D_{1-4}\ (um^2/s)$')
    
    fig, ax = plt.subplots()
    data = all_outputs[all_outputs['sampletype'].str.contains('CTL small', case=False, na=False)]
    p = sns.histplot(x = 'D_1_4', data = data, hue = 'sampletype', element = 'step',
                     stat = 'percent', common_norm = False, \
                         multiple = 'layer', shrink = 1, bins = bins)
    ax.set_xscale('log')
    p.set(xlabel = r'$D_{1-4}\ (um^2/s)$')
    
    fig, ax = plt.subplots()
    data = all_outputs[all_outputs['sampletype'].str.contains('CTL big', case=False, na=False)]
    p = sns.histplot(x = 'D_1_4', data = data, hue = 'sampletype', element = 'step',
                     stat = 'percent', common_norm = False, \
                         multiple = 'layer', shrink = 1, bins = bins)
    ax.set_xscale('log')
    p.set(xlabel = r'$D_{1-4}\ (um^2/s)$')
    
    
    
    
    bins = np.linspace(-0.2, 2.0, num = 23)
    fig, ax = plt.subplots()
    data = all_outputs[all_outputs['sampletype'].str.contains('KO small', case=False, na=False)]
    p = sns.histplot(x = 'alpha', data = data, hue = 'sampletype', element = 'step',
                     stat = 'percent', common_norm = False, \
                         multiple = 'layer', shrink = 1, bins = bins)
    p.set(xlabel = r'$alpha$')

    fig, ax = plt.subplots()
    data = all_outputs[all_outputs['sampletype'].str.contains('KO big', case=False, na=False)]
    p = sns.histplot(x = 'alpha', data = data, hue = 'sampletype', element = 'step',
                     stat = 'percent', common_norm = False, \
                         multiple = 'layer', shrink = 1, bins = bins)
    p.set(xlabel = r'$alpha$')
    
    fig, ax = plt.subplots()
    data = all_outputs[all_outputs['sampletype'].str.contains('CTL small', case=False, na=False)]
    p = sns.histplot(x = 'alpha', data = data, hue = 'sampletype', element = 'step',
                     stat = 'percent', common_norm = False, \
                         multiple = 'layer', shrink = 1, bins = bins)
    p.set(xlabel = r'$alpha$')
    
    fig, ax = plt.subplots()
    data = all_outputs[all_outputs['sampletype'].str.contains('CTL big', case=False, na=False)]
    p = sns.histplot(x = 'alpha', data = data, hue = 'sampletype', element = 'step',
                     stat = 'percent', common_norm = False, \
                         multiple = 'layer', shrink = 1, bins = bins)
    p.set(xlabel = r'$alpha$')

    
    # except Exception as e:
    #     # Handle the error and skip the iteration
    #     print(f"Error processing video {subfolder_name[i]}. Skipping to the next video.")

