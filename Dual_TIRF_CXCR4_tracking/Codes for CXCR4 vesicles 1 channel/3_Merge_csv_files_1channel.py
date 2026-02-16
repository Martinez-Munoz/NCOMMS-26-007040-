# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 10:46:09 2025

@author: rpons
"""

import os
import pandas as pd
import re

# Below, the parameters to check before running the code:
################################################################################################################

# Specify the folder containing the CSV files
folder_path = os.path.normpath(r"/Volumes/CABIMER 5/Spinning disc/20221111 KIN X4")
experiment_to_analyze = "Experiment_1"

################################################################################################################

def get_trajectory_densities(data):
    from itertools import product
    # Extract 'CTL' or 'KO' based on the first appearance in the filename
    def extract_group(filename):
        match = re.search(r'(VCAM|KO)', filename)
        return match.group(0) if match else None
    
    # Map the group (CTL or KO) to each filename
    filename_to_group = {filename: extract_group(filename) for filename in data['filename'].unique()}
    
    # Define the expected combinations for each filename
    expected_combinations = []
    for filename in data['filename'].unique():
        group = filename_to_group[filename]  # Get the group (CTL or KO) for the current filename
        for sampletype in ['big', 'small']:
            expected_combinations.append({
                'filename': filename,
                'sampletype': f"{group} {sampletype}",
            })
    
    expected_combinations_df = pd.DataFrame(expected_combinations)
    
    grouped = data.groupby(['filename', 'sampletype']).size().reset_index(name='count')
    
    # Merge the grouped data with the expected combinations
    result = expected_combinations_df.merge(grouped, on=['filename', 'sampletype'], how='left')
    
    factor_by_filename = data[['filename', 'CELL_AREA']].drop_duplicates()
    result = result.merge(factor_by_filename, on='filename', how='left', suffixes=('', '_from_filename'))
    
    # Fill missing counts with 0
    result['count'] = result['count'].fillna(0)
    
    result['density'] = result['count'] / result['CELL_AREA']
    
    return result

search_path = os.path.join(folder_path, experiment_to_analyze)
file_naming_pattern = "diffusion_parameters_cell"  # The substring to search for in file names

# Initialize an empty list to store individual DataFrames
data_frames = []

# Loop through all files in the experiment folder
for file_name in os.listdir(search_path):
    # Check if the file matches the naming pattern and is a CSV
    if file_naming_pattern in file_name and file_name.endswith(".csv"):
        file_path = os.path.join(search_path, file_name)
        try:
            # Read the CSV file into a DataFrame
            df = pd.read_csv(file_path)
            data_frames.append(df)
        except Exception as e:
            print(f"Error reading {file_name}: {e}")

# Combine all DataFrames into a single DataFrame
if data_frames:
    # Combine the DataFrames directly without introducing NaN
    combined_df = pd.concat(data_frames, ignore_index=True)
    print("CSV files successfully combined!")
else:
    print("No matching CSV files found.")
    combined_df = pd.DataFrame()
    
if not combined_df.empty:
    all_densities = get_trajectory_densities(combined_df)
    
    # Optionally save the combined DataFrame to a new CSV file
    output_path = os.path.join(search_path, f"{experiment_to_analyze}_diffusion_parameters_combined.csv")
    combined_df.to_csv(output_path, index=False)
    output_path = os.path.join(search_path, f"{experiment_to_analyze}_trajectory_densities.csv")
    all_densities.to_csv(output_path, encoding='utf-8', index=False)

    print(f"Combined CSV and densities file saved to {search_path}")
