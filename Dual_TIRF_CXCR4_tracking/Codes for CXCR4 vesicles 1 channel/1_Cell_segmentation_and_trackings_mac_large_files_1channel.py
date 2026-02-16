# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 13:52:37 2024

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
import time
import re
import glob

# Below, the parameters to check before running the code:
################################################################################################################

# In the path below, you should specify where all of your data is stored.
folder_name = os.path.normpath(r"/Volumes/CABIMER 5/Spinning disc/20221111 KIN X4")

# In case you don't want to analyze all videos, below you can specify which experiment, condition, timepoint, and specific video to analyze.
# None --> analyzes all videos
experiment_to_analyze = "Experiment_1"  # Set to "experiment1" or None for all
condition_to_analyze = "VCAM"  # Set to "KO" or None for all
time_point_to_analyze = None  # Set to "10 min" or None for all
video_to_analyze = None       # Set to "video1" or None for all

# Below, you can specify what data you want to export.
cell_mask_export = True #Do you want to export a tif file with the cell segmentation? (True=yes, False=no)
csv_export = True #Do you want to export a csv file with the tracks? (True=yes, False=no)
tif_export = False #Do you want to export a tif file of the original video + overlayed tracks? (True=yes, False=no)

#Specify pixel size
pixel_size = 0.0645

################################################################################################################


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


def import_video(filepath, set_preview, set_roi):
    pattern = os.path.join(filepath, '*.TIF')
    file_list = glob.glob(pattern)
    file_list = sorted(file_list, key=lambda x: int(re.search(r'(\d+)(?=\D*$)', os.path.basename(x)).group()))
    if not set_preview:
        if not np.any(set_roi):
            frames = [imread(file) for file in file_list]
            video = np.asarray(frames)[:,0,:,:]
        if np.any(set_roi):
            first_frame = imread(file_list[0])
            width = first_frame.shape[-1]  # Assume last dimension is width
            frames = [imread(file)[:, set_roi[0,1]:set_roi[1,1], set_roi[0,0]:set_roi[1,0]] for file in file_list]
            video_roi1 = np.asarray(frames)[:,0,:,:]
            frames = [imread(file)[:, set_roi[0,1]:set_roi[1,1], set_roi[0,0]+int(width/2):set_roi[1,0]+int(width/2)] for file in file_list]
            video_roi2 = np.asarray(frames)[:,0,:,:]
            video = np.concatenate((video_roi1, video_roi2), axis=2)
    if set_preview:
        frames = [imread(file) for file in file_list][0]
        video = np.asarray(frames)[0,:,:]
    return video

def adjust_gamma(image, gamma):
	# build a lookup table mapping the pixel values [0, 255] to
	# their adjusted gamma values
    image = image/np.max(image)
    image_gamma = np.power(image, gamma)
    return image_gamma

def cellpose_segmentation(image, total_size):
    # The model chosen is nuclei, as in the images it was performing better than cyto (which was missing detections in images where the cell occupied a large space)
    if total_size < 1 * 1024 * 1024 * 1024:
        model = models.Cellpose(gpu=False, model_type='nuclei')
    else:
        model = models.Cellpose(gpu=False, model_type='cyto')
    cell_mask, flows, styles, diams = model.eval(
        [image],                  # List of images to process
        diameter=200,            # Automatically estimate cell diameter
        channels=[0, 0],          # Specify channels: [cytoplasm, nucleus]
        flow_threshold=0.5,       # Set flow threshold
        cellprob_threshold=0.0    # Set cell probability threshold
    )
    #cell_mask=cell_mask[0][start_y:start_y + image.shape[0], start_x:start_x + image.shape[1]]
    cell_mask = cell_mask[0]
    
    # Cellpose sometimes fails at segmenting cells if these occupy a large area of the whole image. In those cases, the image is embedded into a larger array, and the cellpose segmentation is ran again.
    if np.sum(cell_mask) <= 1:
        big_array_shape = (1000, 1000)
        big_array = np.zeros(big_array_shape, dtype=image.dtype)+image[0,0]
        
        # Calculate the start indices for centering the smaller array
        start_y = (big_array_shape[0] - image.shape[0]) // 2
        start_x = (big_array_shape[1] - image.shape[1]) // 2
        
        # Insert the smaller array into the larger array
        big_array[start_y:start_y + image.shape[0], start_x:start_x + image.shape[1]] = image
        
        # plt.figure()
        # plt.imshow(big_array)
        
        model = models.Cellpose(gpu=True, model_type='nuclei')
        cell_mask, flows, styles, diams = model.eval(
            [big_array],                  # List of images to process
            diameter=200,            # Automatically estimate cell diameter
            channels=[0, 0],          # Specify channels: [cytoplasm, nucleus]
            flow_threshold=0.5,       # Set flow threshold
            cellprob_threshold=0.0    # Set cell probability threshold
        )
        #cell_mask=cell_mask[0][start_y:start_y + image.shape[0], start_x:start_x + image.shape[1]]
        cell_mask = cell_mask[0][start_y:start_y + image.shape[0], start_x:start_x + image.shape[1]]
        
    
    # Initialize a label array to store watershed results
    final_labels = np.zeros_like(cell_mask, dtype=np.int32)
    
    # Loop through each unique cell ID in the `cell_mask`
    unique_cells = np.unique(cell_mask)
    for cell_id in unique_cells:
        if cell_id == 0:  # Skip the background
            continue
        
        # Extract the binary mask for the current cell
        single_cell_mask = (cell_mask == cell_id).astype(np.uint8)
        
        # Compute the distance transform for the single cell mask
        distance = ndi.distance_transform_edt(single_cell_mask)
        distance_noisy = distance + np.random.uniform(low=-1e-6, high=1e-6, size=distance.shape)
        
        # Find local maxima for markers
        coords = peak_local_max(distance_noisy, footprint=np.ones((250, 250)), labels=single_cell_mask)
        marker_mask = np.zeros(distance.shape, dtype=bool)
        marker_mask[tuple(coords.T)] = True
        markers, _ = ndi.label(marker_mask)
        
        # Apply Watershed Segmentation on the single cell mask
        cell_labels = watershed(-distance, markers, mask=single_cell_mask)
        
        # Add the segmented cell to the final labels
        final_labels[cell_labels > 0] = cell_id * 1000 + cell_labels[cell_labels > 0]
    
    # plt.figure()
    # plt.imshow(flows[0][2]) 
    # plt.colorbar(label='Cell Probability')
    
    # plt.figure()
    # plt.imshow(final_labels)
    
    return final_labels

def retrieve_intensities(video, frame, localizations):
    mean_intensities = np.empty((0,1))
    median_intensities = np.empty((0,1))
    max_intensities = np.empty((0,1))
    min_intensities = np.empty((0,1))
    total_intensities = np.empty((0,1))
    std_intensities = np.empty((0,1))
    for i in range(len(localizations)):
        radius = int(localizations[i,2])  # Approximate radius of the blob
    
        # Extract a square region around the blob
        min_row = max(int(localizations[i,0]) - radius, 0)
        max_row = min(int(localizations[i,0]) + radius, video[i].shape[0])
        min_col = max(int(localizations[i,1]) - radius, 0)
        max_col = min(int(localizations[i,1]) + radius, video[i].shape[1])      
        # Extract the region and calculate the integrated intensity
        region = video[i][min_row:max_row, min_col:max_col]
        
        mean_intensities = np.vstack([mean_intensities, np.mean(region)])
        median_intensities = np.vstack([median_intensities, np.median(region)])
        min_intensities = np.vstack([min_intensities, np.min(region)])
        max_intensities = np.vstack([max_intensities, np.max(region)])
        total_intensities = np.vstack([total_intensities, np.sum(region)])
        std_intensities = np.vstack([std_intensities, np.std(region)])
    
    return mean_intensities, median_intensities, min_intensities, max_intensities, total_intensities, std_intensities

def detect_blobs(video, min_sigma, max_sigma, threshold):
    localizations_blobs = []
    for frame in range(video.shape[0]):
        #print(f"Detecting blobs in frame {frame}...")
        blobs_log = feature.blob_log(video[frame], min_sigma=min_sigma,  max_sigma=max_sigma, threshold=threshold)
        #blobs_log = feature.blob_doh(video[frame], min_sigma=min_sigma,  max_sigma=max_sigma, threshold=0.005)
        mean_intensities, median_intensities, min_intensities, max_intensities, total_intensities, std_intensities = retrieve_intensities(video, frame, blobs_log)
        localizations_blobs.append(np.hstack((blobs_log, mean_intensities, median_intensities, min_intensities, max_intensities, total_intensities, std_intensities)))
    return localizations_blobs

def detect_rings(video, dp, minDist, param1, param2, minRadius, maxRadius):
    #laplacian = cv2.Laplacian(video, cv2.CV_64F, ksize=9)
    laplacian = np.empty_like(video).astype(np.uint8)
    # Process each frame
    for i in range(video.shape[0]):  # Iterate over frames
        frame = video[i]  # Extract the ith frame (2D grayscale)
        # Apply Laplacian filter
        laplacian_frame = cv2.Laplacian(frame, cv2.CV_64F, ksize=9)
        laplacian_frame = np.max(laplacian_frame, keepdims=True)-laplacian_frame
        laplacian_frame = laplacian_frame/np.max(laplacian_frame, keepdims=True)
        laplacian_frame = ((laplacian_frame) * 255).astype(np.uint8)
        laplacian_frame = cv2.GaussianBlur(laplacian_frame, ksize=(9, 9), sigmaX=3)
        # Store the processed frame
        laplacian[i] = laplacian_frame
    
    localizations_rings = []
    for i in range(video.shape[0]):
        rings = cv2.HoughCircles(
        laplacian[i], 
        cv2.HOUGH_GRADIENT_ALT, 
        dp=dp, 
        minDist=minDist, 
        param1=param1, 
        param2=param2, 
        minRadius=minRadius, 
        maxRadius=maxRadius
        )
        
        # circles = rings.copy() if rings is not None else None
        
        if rings is None:
            localizations_rings.append(np.empty((0, 9), dtype=np.float32))
        else:
            rings = rings[0]
            rings[:, [0, 1]] = rings[:, [1, 0]]
            
            mean_intensities, median_intensities, min_intensities, max_intensities, total_intensities, std_intensities = retrieve_intensities(video, i, rings)

            localizations_rings.append(np.hstack((rings, mean_intensities, median_intensities, min_intensities, max_intensities, total_intensities, std_intensities)))
        
        # if circles is not None:
        #     # Convert circles to integers
        #     circles = np.uint16(np.around(circles))
        
        #     # Copy the original frame to draw on it
        #     output_image = cv2.cvtColor(laplacian[i], cv2.COLOR_GRAY2BGR)
        
        #     # Draw each detected circle
        #     for circle in circles[0, :]:
        #         x, y, radius = circle
        #         # Draw the circle's perimeter
        #         cv2.circle(output_image, (x, y), radius, (0, 255, 0), 2)
        #         # Draw the circle's center
        #         cv2.circle(output_image, (x, y), 2, (0, 0, 255), 3)
        
        #     # Display the resulting image with circles
        #     plt.figure(figsize=(10, 10))
        #     plt.imshow(cv2.cvtColor(output_image, cv2.COLOR_BGR2RGB))
        #     plt.title(f"Detected Circles {i}")
        #     plt.axis("off")
        #     plt.show()
        # else:
        #     print("No circles were detected.") 
    return localizations_rings

def ensure_rings(trajectories, video, min_length, threshold):
    for particle in trajectories.particle.unique():
        ring_traj = trajectories[trajectories["particle"] == particle]
        if len(ring_traj) < min_length:
            trajectories = trajectories[trajectories["particle"] != particle]
        else:
            r=int(np.mean(ring_traj.radius))
            inset = []
            inset_inner = []
            for frame in ring_traj.frame:
                #print(frame)
                xr = int(ring_traj.x[ring_traj["frame"] == frame].iloc[0])
                yr = int(ring_traj.y[ring_traj["frame"] == frame].iloc[0])
                inset_frame = video[frame, xr-r:xr+r, yr-r:yr+r]
                inset_inner_frame = video[frame, xr-2:xr+2, yr-2:yr+2]
                #plt.figure()
                #plt.imshow(inset_frame)
                inset.append(inset_frame)
                inset_inner.append(inset_inner_frame)
            inset = np.array(inset)
            inset_average = np.mean(inset, axis=0)
            inset_inner = np.array(inset_inner)
            inset_inner_average = np.mean(inset_inner, axis=0)
            #plt.figure()
            #plt.imshow(inset_average)
            if np.mean(inset_inner_average) >= np.percentile(inset_average, threshold):
                trajectories = trajectories[trajectories["particle"] != particle]
    return trajectories

def ring_corrections(trajectories_blobs, trajectories_rings, addRings):
    
    def majority_positive_or_negative(group):
        # Count positive values in the group
        positive_count = (group > 0).sum()
        # Count negative values in the group
        negative_count = (group < 0).sum()
        
        # Compare the counts
        if positive_count > negative_count:
            return "Positive"
        elif negative_count > positive_count:
            return "Negative"
        else:
            return "Equal"
    
    temp = trajectories_rings.copy()
    #temp["candidates particle"] =  ""
    temp["candidate particles"] = [[] for _ in range(len(temp))]
    #temp=[]
    candidates_exclusion = []
    #for i in range(len(trajectories_rings)):
    for idx in trajectories_rings.index:
        xr = np.array(trajectories_rings.loc[idx, trajectories_rings.columns[0]])
        yr = np.array(trajectories_rings.loc[idx, trajectories_rings.columns[1]])
        r = np.array(trajectories_rings.loc[idx, trajectories_rings.columns[2]])
        frame = np.array(trajectories_rings.loc[idx, trajectories_rings.columns[3]])
        candidates = trajectories_blobs[trajectories_blobs["frame"] == frame].copy()
        candidates_coord = np.array(candidates.iloc[:,0:2])
        distances = np.sqrt(np.sum((candidates_coord - np.array([[xr, yr]])) ** 2, axis=1))
        candidates["distance"] = distances
        candidates["distance vs radii"] = distances-r/2
        candidates = candidates[candidates["distance"] <= r+2]
        #temp.iloc[i, 5] = np.array2string(candidates["particle"].to_numpy())
        #temp.iloc[i, 5] = [candidates["particle"].to_numpy()]
        temp.at[idx, "candidate particles"] = candidates["particle"].to_numpy()
        candidates_exclusion.append(candidates)
    if candidates_exclusion:
        candidates_exclusion = pd.concat(candidates_exclusion, ignore_index=True)
    
        grouped = candidates_exclusion.groupby("particle")["distance vs radii"]
        majority_result = grouped.apply(majority_positive_or_negative)
        positive_particles = majority_result[majority_result == "Positive"].index.tolist()
        trajectories_blobs = trajectories_blobs[~trajectories_blobs["particle"].isin(positive_particles)]
        
        if addRings:
            for particle in temp.particle.unique():
                candidates_particles = temp[temp["particle"] == particle]["candidate particles"]
                candidates_particles = np.concatenate([np.array(item).ravel() for item in candidates_particles])
                if np.all(np.isin(candidates_particles, positive_particles)):
                    ring_traj = trajectories_rings[trajectories_rings["particle"] == particle].copy()
                    ring_traj["particle"] = np.max(trajectories_blobs["particle"])+1
                    trajectories_blobs = pd.concat([trajectories_blobs, ring_traj], ignore_index=True)
                
    return trajectories_blobs



def linking(localizations, search_range, allow_n_gaps):
    logging.getLogger('trackpy').setLevel(logging.WARNING)
    localizations_dataframe = []
    # Loop through each frame's data and add a frame column
    for frame_number, coords in enumerate(localizations):
        # Create a DataFrame for each frame's localizations
        df_frame = pd.DataFrame(coords[:,0:3], columns=['x', 'y', 'radius'])
        df_frame['radius'] = ((df_frame['radius']-2)/0.55555556)*10**-10
        df_frame['frame'] = frame_number  # Add frame number
        localizations_dataframe.append(df_frame)
    # Concatenate all frames into a single DataFrame
    df_all_frames = pd.concat(localizations_dataframe, ignore_index=True)
    search_range[0]=search_range[0]*10**-10
    trajectories = tp.link_df(df_all_frames, search_range=search_range, memory=1, pos_columns=['radius','y','x'])
    trajectories['radius'] = (trajectories['radius']/10**-10)*0.55555556+2
    return trajectories


def classify(trajectories, radius_threshold):
    avg_radius_per_particle = (trajectories.groupby('particle')['radius'].mean())*pixel_size*np.sqrt(2)
    
    avg_radius_df = avg_radius_per_particle.reset_index()
    
    # Assign a size category based on the average radius
    avg_radius_df['size_category'] = avg_radius_df['radius'].apply(
        lambda x: 'small' if x < radius_threshold else 'big'
    )
    
    # Merge the size category back into the original trajectories DataFrame
    trajectories_with_size = trajectories.merge(
        avg_radius_df[['particle', 'size_category']], on='particle'
    )
    
    # Split the DataFrame into small and big based on size_category
    trajectories_small = trajectories_with_size[
        trajectories_with_size['size_category'] == 'small'
    ].drop(columns=['size_category'])
    
    trajectories_big = trajectories_with_size[
        trajectories_with_size['size_category'] == 'big'
    ].drop(columns=['size_category'])
    
    return trajectories_small, trajectories_big
   
def writeheader(csvwriter):
    # Define the header rows with Unicode for 'µ'
    line1 = [
        'LABEL', 'ID', 'TRACK_ID', 'QUALITY', 'POSITION_X', 'POSITION_Y',
        'POSITION_Z', 'POSITION_T', 'FRAME', 'RADIUS', 'VISIBILITY',
        'MANUAL_SPOT_COLOR', 'MEAN_INTENSITY_CH1', 'MEDIAN_INTENSITY_CH1',
        'MIN_INTENSITY_CH1', 'MAX_INTENSITY_CH1', 'TOTAL_INTENSITY_CH1',
        'STD_INTENSITY_CH1', 'CONTRAST_CH1', 'SNR_CH1'
    ]
    line2 = [
        'Label', 'Spot ID', 'Track ID', 'Quality', 'X', 'Y', 'Z', 'T',
        'Frame', 'Radius', 'Visibility', 'Manual spot color', 'Mean intensity ch1',
        'Median intensity ch1', 'Min intensity ch1', 'Max intensity ch1',
        'Sum intensity ch1', 'Std intensity ch1', 'Contrast ch1',
        'Signal/Noise ratio ch1'
    ]
    line3 = [
        'Label', 'Spot ID', 'Track ID', 'Quality', 'X', 'Y', 'Z', 'T',
        'Frame', 'R', 'Visibility', 'Spot color', 'Mean ch1', 'Median ch1',
        'Min ch1', 'Max ch1', 'Sum ch1', 'Std ch1', 'Ctrst ch1', 'SNR ch1'
    ]
    line4 = [
        '', '', '', '(quality)', '(\u03bcm)', '(\u03bcm)', '(\u03bcm)', '(frame)', '',
        '(\u03bcm)', '', '', '(counts)', '(counts)', '(counts)', '(counts)', '(counts)',
        '(counts)', '', ''
    ]
    
    #line1=[x.encode('cp1252') for x in line1]
    #line2=[x.encode('cp1252') for x in line2]
    #line3=[x.encode('cp1252') for x in line3]
    #line4=[x.encode('cp1252') for x in line4]
    
    # Write header rows
    csvwriter.writerow(line1)
    csvwriter.writerow(line2)
    csvwriter.writerow(line3)
    csvwriter.writerow(line4)

def export_tracks(trajectories, localizations, csv_track_path):
    #print('Tracks file:', csv_track_path)
    
    # Open the CSV file in write mode with UTF-8 encoding
    with open(csv_track_path, 'w', encoding='utf-8-sig', newline='') as csv_track_file:
        spotswriter = csv.writer(csv_track_file)
        
        # Assuming writeheader is a defined function that writes the CSV header
        writeheader(spotswriter)
        
        # Iterate through unique particle IDs in the trajectories DataFrame
        for track_id, trajectory in trajectories.groupby('particle'):
            indices = trajectory.index
            #print(trajectory)
            #print(particle_id)
            points = trajectory[['x', 'y', 'frame']].values
            if len(points) >= 10:
                for j in range(len(points)):
                    #print(j)
                    spot_id = indices[j]
                    x=points[j,1]
                    y=points[j,0]
                    frame=points[j,2]
                    q=np.nan
                    snr=np.nan
                    mean=localizations[int(frame)][(localizations[int(frame)][:, 0] == y) & (localizations[int(frame)][:, 1] == x)][0,3]
                    median=localizations[int(frame)][(localizations[int(frame)][:, 0] == y) & (localizations[int(frame)][:, 1] == x)][0,4]
                    mmin=localizations[int(frame)][(localizations[int(frame)][:, 0] == y) & (localizations[int(frame)][:, 1] == x)][0,5]
                    mmax=localizations[int(frame)][(localizations[int(frame)][:, 0] == y) & (localizations[int(frame)][:, 1] == x)][0,6]
                    total=localizations[int(frame)][(localizations[int(frame)][:, 0] == y) & (localizations[int(frame)][:, 1] == x)][0,7]
                    std=localizations[int(frame)][(localizations[int(frame)][:, 0] == y) & (localizations[int(frame)][:, 1] == x)][0,8]
                    radius=localizations[int(frame)][(localizations[int(frame)][:, 0] == y) & (localizations[int(frame)][:, 1] == x)][0,2]
                    #print spot.getFeatures() # show all spot features
                    spotlist = ['ID'+str(spot_id), spot_id, track_id, 0, x*pixel_size, \
                    y*pixel_size, 0, frame, \
                    frame, radius*pixel_size, 1,\
                    'None', mean, median, \
                    mmin, mmax, \
                    total, std, 0, \
                    0]
                    
                    #print spotlist
                    #model.getLogger().log('\tspot ID = ' + str(sid) + ': x='+str(x)+', y='+str(y)+', t='+str(t)+', q='+str(q) + ', snr='+str(snr) + ', mean = ' + str(mean))
                    #spotlist=[str(x).encode('cp1252') for x in spotlist]
                    spotswriter.writerow(spotlist)     
    csv_track_file.close()    
    
def export_video(trajectories, video, tif_video_path):
    plt.rcParams['figure.dpi'] = 600

    particle_ids = trajectories['particle'].unique()
    color_map = {
        particle_id: (random.randint(0, 255), random.randint(0, 255), random.randint(0, 255))
        for particle_id in particle_ids
    }

    fade_duration = 15  # Number of frames to fade out

    multi_channel_frames = []

    for frame_number, frame in enumerate(video):
        # Scale the original 16-bit frame to 8-bit
        frame_8bit = (frame/np.max(video)*255).astype(np.uint8)

        # Create a black-and-white RGB frame by copying grayscale data to all channels
        bw_frame = np.stack([frame_8bit] * 3, axis=-1)

        # Create a transparent overlay for drawing the fading trajectories
        overlay = np.zeros((*bw_frame.shape[:2], 4), dtype=np.uint8)  # 4-channel RGBA image

        # Draw trajectories with fading effect
        particles_in_frame = trajectories[trajectories['frame'] <= frame_number]
        for particle_id, trajectory in particles_in_frame.groupby('particle'):
            color = color_map[particle_id]  # Get base color for this particle
            points = trajectory[['x', 'y', 'frame']].values
            if len(points) >= 12:
                for j in range(1, len(points)):
                    x0, y0, frame0 = int(points[j - 1][0]), int(points[j - 1][1]), int(points[j - 1][2])
                    x1, y1, frame1 = int(points[j][0]), int(points[j][1]), int(points[j][2])
    
                    # Calculate fade factor based on how many frames back the segment is
                    frames_ago = frame_number - frame1
                    fade_factor = max(0, 1 - frames_ago / fade_duration)  # Linear fade out to 0
    
                    # Apply fade factor to the color and set alpha
                    faded_color = tuple(int(c) for c in color) + (int(255 * fade_factor),)
                    
                    # Draw the line with the faded color on the transparent overlay
                    cv2.line(overlay, (y0, x0), (y1, x1), faded_color, thickness=1, lineType=cv2.LINE_AA)

        # Blend the overlay with the base black-and-white frame
        overlay_pil = Image.fromarray(overlay, mode="RGBA")
        bw_frame_pil = Image.fromarray(bw_frame, mode="RGB")
        blended_image = Image.alpha_composite(bw_frame_pil.convert("RGBA"), overlay_pil)

        # Convert back to RGB for saving
        final_image = blended_image.convert("RGB")
        
        # #Draw particle IDs on the final image
        # draw = ImageDraw.Draw(final_image)
        # # Get particles in the current frame
        # current_particles = trajectories[trajectories['frame'] == frame_number]
        # for idx, row in current_particles.iterrows():
        #     x = int(row['x'])
        #     y = int(row['y'])
        #     particle_id = row['particle']
        #     color = color_map[particle_id]
        #     text = str(particle_id)
        #     text_position = (y + 5, x + 5)  # Adjust as needed
        #     draw.text(text_position, text, fill=tuple(color))
        
        
        multi_channel_frames.append(final_image)

    # Save as a multi-page TIFF
    multi_channel_frames[0].save(tif_video_path, save_all=True, append_images=multi_channel_frames[1:], compression="tiff_deflate")


def analyze_video(video, total_size):
    #Splitting the video into the two channels
    video_cxcr4 = video
    #video_cxcr4 = video_cxcr4[0:100]
    #video_rab = video_rab[0:100]
    
    #Segmenting cells using Cellpose on the time-averaged Rab channel
    image_cells = adjust_gamma(np.mean(video_cxcr4, axis=0), 1) #a gamma filter is applied to help segment cells of varying expression levels
    # plt.figure()
    # plt.imshow(image_cells)
    segmented_cells = cellpose_segmentation(image_cells, total_size)
    
    #Detecting and tracking CXCR4 and Rab structures. Note: the signal in all frames is normalized to 1 in order to (1) keep the detections consistent across cells of varying expression levels, and (2) correct for photobleaching.
    video_cxcr4_norm = video_cxcr4/np.max(video_cxcr4, axis=(1,2), keepdims=True)
    localizations_blobs_cxcr4 = detect_blobs(video_cxcr4_norm, min_sigma=2, max_sigma=8, threshold=0.04000)
    localizations_rings_cxcr4 = detect_rings(video_cxcr4_norm, dp=1.4, minDist=8, param1=1, param2=0.6, minRadius=2, maxRadius=15)
    
    return segmented_cells, localizations_blobs_cxcr4, localizations_rings_cxcr4


def linking_general(localizations_blobs_cxcr4, localizations_rings_cxcr4, video):
    #Splitting the video into the two channels
    video_cxcr4 = video
    video_cxcr4_norm = video_cxcr4/np.max(video_cxcr4, axis=(1,2), keepdims=True)
    
    trajectories_blobs_cxcr4 = linking(localizations_blobs_cxcr4, search_range=[3,5,5], allow_n_gaps=1) #search range determines maximum frame-to-frame changes in the radius (Ar), y (Ay), and x (Ax) --> [Ar, Ay, Ax]. All units in pixels.
    trajectories_rings_cxcr4 = linking(localizations_rings_cxcr4, search_range=[3,5,5], allow_n_gaps=3)
    trajectories_rings_cxcr4 = ensure_rings(trajectories_rings_cxcr4, video_cxcr4_norm, min_length=3, threshold=80)
    trajectories_cxcr4 = ring_corrections(trajectories_blobs_cxcr4, trajectories_rings_cxcr4, addRings=False)
    
    return trajectories_cxcr4
    

def correct_coordinates(l1, l2, df1, roi):
    constants_per_list = [[roi[0,1], roi[0,0]], [roi[0,1], roi[0,0]]]
    # Combine the lists and their corresponding constants
    lists = [l1, l2]
    # Add the corresponding constants to each column in the arrays
    for array_list, constants in zip(lists, constants_per_list):
        for arr in array_list:
            # Add constants to each column
            arr[:, :len(constants)] += constants  # Only modify columns up to the length of constants
            
    constants_list = [{'x': roi[0,1], 'y': roi[0,0]}, {'x': roi[0,1], 'y': roi[0,0]}]
    dataframes = [df1]
    for df, constants in zip(dataframes, constants_list):
        for col, constant in constants.items():
            df[col] += constant  # Add the constant to the specified column
        
    return l1, l2, df1


def ensure_unique_ids(list_dfs):
    current_max_id = 0
    adjusted_dataframes = []
    for df in list_dfs:
        df['particle'] += current_max_id+1  # Adjust IDs by adding the current max
        current_max_id = df['particle'].max()  # Update the current max ID
        adjusted_dataframes.append(df)
    return list_dfs




for i in range(len(video_folder_paths)):
    relative_path = os.path.relpath(video_folder_paths[i], folder_name)
    relative_name = relative_path.split(os.sep)[0] + "_" + relative_path.split(os.sep)[1] + "_" + relative_path.split(os.sep)[2] + "_" + relative_path.split(os.sep)[3]
    start_time = time.perf_counter()
    print(f"({i+1}/{len(video_folder_paths)}) Analyzing video {relative_path}...")
    try:
        #Importing video
        filepath = video_folder_paths[i]
        total_size = sum(os.path.getsize(f) for f in glob.glob(os.path.join(filepath, '*.TIF')))
        pattern = os.path.join(filepath, '*.TIF')
        file_list = glob.glob(pattern)
        file_list = sorted(file_list, key=lambda x: int(re.search(r'(\d+)(?=\D*$)', os.path.basename(x)).group()))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            if total_size < 1 * 1024 * 1024 * 1024:  # 1 GB in bytes
                video = import_video(filepath, set_preview=False, set_roi=False)
                segmented_cells, localizations_blobs_cxcr4, localizations_rings_cxcr4 = analyze_video(video, total_size)
                trajectories_cxcr4 = linking_general(localizations_blobs_cxcr4, localizations_rings_cxcr4, video)
            else:
                print("This is a large file (>1GB). Automatic ROIs will be used.")
                preview = import_video(filepath, set_preview=True, set_roi=False)
                preview_image_cells = adjust_gamma(preview, 0.5)
                segmented_cells = cellpose_segmentation(preview_image_cells, total_size)
                unique_values, counts = np.unique(segmented_cells, return_counts=True)
                values_to_keep = unique_values[counts >= 1000]
                localizations_blobs_cxcr4 = []
                localizations_rings_cxcr4 = []
                trajectories_cxcr4 = []
                for i in values_to_keep[1:]:
                    nonzero_indices = np.nonzero(segmented_cells == i)
                    y_min, y_max = np.min(nonzero_indices[0]), np.max(nonzero_indices[0])
                    x_min, x_max = np.min(nonzero_indices[1]), np.max(nonzero_indices[1])
                    roi_lims = np.array([[x_min, y_min], [x_max, y_max]])
                    video_roi = import_video(filepath, set_preview=False, set_roi=roi_lims)
                    segmented_cells_roi, localizations_blobs_cxcr4_roi, localizations_rings_cxcr4_roi = analyze_video(video_roi, total_size)
                    trajectories_cxcr4_roi = linking_general(localizations_blobs_cxcr4_roi, localizations_rings_cxcr4_roi, video_roi)
                    localizations_blobs_cxcr4_roi, localizations_rings_cxcr4_roi, trajectories_cxcr4_roi = correct_coordinates(localizations_blobs_cxcr4_roi, localizations_rings_cxcr4_roi, trajectories_cxcr4_roi, roi_lims)
                    localizations_blobs_cxcr4.append(localizations_blobs_cxcr4_roi)
                    localizations_rings_cxcr4.append(localizations_rings_cxcr4_roi)
                    trajectories_cxcr4.append(trajectories_cxcr4_roi)
                localizations_blobs_cxcr4 = [np.vstack(arrays) for arrays in zip(*localizations_blobs_cxcr4)]
                localizations_rings_cxcr4 = [np.vstack(arrays) for arrays in zip(*localizations_rings_cxcr4)]
                trajectories_cxcr4 = ensure_unique_ids(trajectories_cxcr4)
                trajectories_cxcr4 = pd.concat(trajectories_cxcr4, axis=0, ignore_index=True)
        
            #Exporting the analysis
            relative_experiment = relative_path.split(os.sep)[0]
            relative_name = relative_path.split(os.sep)[0] + "_" + relative_path.split(os.sep)[1] + "_" + relative_path.split(os.sep)[2] + "_" + relative_path.split(os.sep)[3]
            save_path = os.path.join(folder_name, relative_experiment, relative_name)
            if cell_mask_export:
                tif_ending = '_Cellpose_cell_mask.tif'
                tif_video_path = save_path + tif_ending
                cell_mask = Image.fromarray(segmented_cells)
                cell_mask.save(tif_video_path, compression="tiff_deflate")
            
            if csv_export:
                csv_ending = '_Tracks_CXCR4.csv'
                csv_track_path = save_path + csv_ending
                # The for loop below converts the sigma of blob_log to the radius by multiplying it by a factor of 2^(1/2).
                for arr in localizations_blobs_cxcr4:
                    arr[:, 2] *= np.sqrt(2)
                localizations_all_cxcr4 = [np.vstack((a, b)) for a, b in zip(localizations_blobs_cxcr4, localizations_rings_cxcr4)]
                export_tracks(trajectories_cxcr4, localizations_all_cxcr4, csv_track_path)
        
            if total_size < 1 * 1024 * 1024 * 1024:
                if tif_export:
                    tif_ending = '_video_CXCR4_overlayed_trajectories.tif'
                    tif_video_path = save_path + tif_ending
                    export_video(trajectories_cxcr4, video_cxcr4[:], tif_video_path)
                    
                    tif_ending = '_video_Rab_overlayed_trajectories.tif'
                    tif_video_path = save_path + tif_ending
                    export_video(trajectories_rab, video_rab[:], tif_video_path)
                
    except Exception as e:
        # Handle the error and skip the iteration
        print(f"Error processing video {relative_name}. Skipping to the next video.")

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time} seconds")
        