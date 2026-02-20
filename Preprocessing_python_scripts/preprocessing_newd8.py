#!/usr/bin/env python3

from pyproj import Transformer
import pandas as pd
from funmixer import (
    check_d8,
    get_sample_graph,
    plot_network,
    set_d8_boundaries_to_zero,
    snap_to_drainage,
)
import matplotlib.pyplot as plt
import numpy as np
import os
import rasterio


"""
Step 1:
This part of the script merges two excel sheets
"""

# Load the two CSV files into dataframes
phosphate_df = pd.read_csv("/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Original/Phosphate_date.csv")
sites_df = pd.read_csv("/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Original/Sites.csv")

# Merge the dataframes on 'Site ID' column
merged_df = pd.merge(phosphate_df, sites_df, on='Site ID', how='left')

# Reorder columns to match the desired order
merged_df = merged_df[[
    'Sample ID', 'Date', 'Time', 'Phosphate (Hanna)', 'Site ID', 'Phosphate st', 'Hanna LR phos', 
    'Org', 'Name in EC', 'Latitude', 'Longitude', 'Watercourse (if known)', 'Catchment (or "tbc")'
]]

# Save the merged dataframe to a new CSV file
merged_df.to_csv("/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Merged_Phosphate_Sites.csv", index=False)

print("Merge complete. Saved as 'Merged_Phosphate_Sites.csv'")

"""
Step 2:
This part of the script adds x and y coordinates (British National Grid) to the merged CSV
"""

# Step 1: Read the merged CSV
merged_df = pd.read_csv("/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Merged_Phosphate_Sites.csv")

# Step 2: Ensure column names are clean
merged_df.columns = merged_df.columns.str.strip()

# Step 3: Define transformer for WGS84 -> British National Grid
transformer = Transformer.from_crs("EPSG:4326", "EPSG:27700", always_xy=True)

# Step 4: Apply transformation to all rows
def latlon_to_xy(lon, lat):
    x, y = transformer.transform(lon, lat)
    return x, y

merged_df[['x', 'y']] = merged_df.apply(
    lambda row: latlon_to_xy(row['Longitude'], row['Latitude']),  
    axis=1,
    result_type='expand'
)

# Step 5: Save the new CSV
merged_df.to_csv("/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Merged_Phosphate_Sites_XY.csv", index=False)

print("✅ Added x and y coordinates (British National Grid) and saved as 'Merged_Phosphate_Sites_XY.csv'")

"""
Step 3:
Filter sites in merged CSV based on date, valid coordinates, non-empty phosphate values, and minimum sample count per site.
"""

# Step 1: Read the merged CSV
merged_df = pd.read_csv("/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Merged_Phosphate_Sites_XY.csv")

# Step 2: Clean column names
merged_df.columns = merged_df.columns.str.strip()

# Step 3: Convert 'Date' column to datetime format
merged_df['Date'] = pd.to_datetime(merged_df['Date'], errors='coerce')

# Step 4: Filter by date range (last two years: 2023-10-19 to 2025-10-19)
start_date = pd.Timestamp('2023-10-19')
end_date = pd.Timestamp('2025-10-19')
filtered_df = merged_df[(merged_df['Date'] >= start_date) & (merged_df['Date'] <= end_date)]

# Step 5: Remove rows with missing or invalid Latitude/Longitude values
# Convert to numeric to handle cases where they're stored as strings
filtered_df['Latitude'] = pd.to_numeric(filtered_df['Latitude'], errors='coerce')
filtered_df['Longitude'] = pd.to_numeric(filtered_df['Longitude'], errors='coerce')

# Step 6: Drop rows where Latitude or Longitude are missing (NaN)
filtered_df = filtered_df.dropna(subset=['Latitude', 'Longitude'])

# Step 7: Remove rows where 'Hanna LR phos' is empty (NaN)
filtered_df = filtered_df.dropna(subset=['Hanna LR phos'])

# Step 8: Filter Site IDs with at least 50 samples
site_counts = filtered_df['Site ID'].value_counts()
sites_with_50_samples = site_counts[site_counts >= 50].index
filtered_df = filtered_df[filtered_df['Site ID'].isin(sites_with_50_samples)]

# Step 9: Save the filtered dataframe to a new CSV
output_path = "/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Filtered_Merged_Phosphate_Sites_XY.csv"
filtered_df.to_csv(output_path, index=False)

# Step 10: Print summary info
unique_sites = filtered_df['Site ID'].nunique()

print(f"✅ Filtered data saved as:\n{output_path}")
print(f"📊 Rows after filtering: {len(filtered_df)}")
print(f"📍 Unique sites after filtering: {unique_sites}")

"""
Step 4:
This part of the script creates a new CSV with unique Site IDs and their x, y coordinates
"""
# Step 1: Load your original CSV
input_csv = "/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Filtered_Merged_Phosphate_Sites_XY.csv"
df = pd.read_csv(input_csv)

# Step 2: Select only the desired columns
new_df = df[['Site ID', 'x', 'y']]

# Step 3: ensure Sample ID is string
new_df['Site ID'] = new_df['Site ID'].astype(str)

# Step 4: Drop duplicates based on Sample ID, keeping the first occurrence
unique_df = new_df.drop_duplicates(subset='Site ID', keep='first')

# Step 5: Save to a new CSV
output_csv = "/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Unique_SiteID_XY.csv"
unique_df.to_csv(output_csv, index=False)

print(f"New CSV saved as {output_csv} with {len(unique_df)} unique rows.")



"""
Step 5
This script demonstrates the preprocessing capabilities of funmixer. Specifically: 
1. Checking that a D8 flow direction raster is correctly formatted for use in funmixer
2. Fixing a D8 raster that has incorrect boundary conditions 
3. Snapping sample sites to the nearest drainage network
"""

### Checking D8 flow directions and fixing boundary conditions ###
# The check_d8 function checks two things.
# 1. if the flow-direction values are the expected values e.g., 0, 1, 2, 4, 8, 16, 32, 64, 128
# 2. if the boundary conditions are correct, i.e., all boundary cells are sinks (0's)
# The example file "d8_bad_bounds.tif" has correct values but incorrect boundary conditions.
# We can test this using the check_d8 function.

check_d8("/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Original/welsh_d8.nc")

# We can fix the boundary conditions using the set_d8_boundaries_to_zero function which sets all boundary cells to 0,
# writing the corrected raster to a new file.
set_d8_boundaries_to_zero("/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Original/welsh_d8.nc")

# Now we can check the corrected raster.
check_d8("welsh_d8_fix_bounds.tif")

### Snapping misaligned sample sites to the nearest drainage network ###
# In general, sample sites are not perfectly aligned with the drainage network, due to
# uncertain locations or the inherent simplification in representing flow using D8. This
# means that when genreating a "sample_network" using the get_sample_graphs function,
# the generated network may be incorrect. e.g., the sample sites may be connected to the wrong
# tributary or simply be disconnected from the network entirely. This can be fixed by snapping
# the sample sites to the nearest drainage network using the snap_to_drainage function.

# Load in real samples
samples = pd.read_csv("/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Unique_SiteID_XY.csv")
print(samples.columns)
sample_x, sample_y = samples["x"], samples["y"]

# Load sample network
sample_network, labels = get_sample_graph(
    flowdirs_filename="welsh_d8_fix_bounds.tif",
    sample_data_filename="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Unique_SiteID_XY.csv",
)

plt.title("Disconnected Network Due to Misaligned Samples")
plot_network(sample_network)
plt.show()

# We can snap the noisy samples to the nearest drainage network using the snap_to_drainage function.
# This requires the flow direction raster and sample sites as input. The drainage_area_threshold parameter
# tells the function to only snap to drainage pixels with drainage area greater than the threshold. This is
# in the same units as the flow direction raster. The plot and save parameters control whether the snapped
# sample sites are plotted and saved to file respectively.
# The nudges parameter allows for manual nudges of sample sites to the nearest drainage pixel. For example, the
# below code nudges the sample site "CG001" by 1000m in the x-direction and -1000m in the y-direction. These
# can be visualised by setting plot=True. It may take some trial and error to get the nudges right and the samples
# snapped to the correct drainage pixel.

print(samples["Site ID"].unique()) 


scale = 10  # metres

nudges = {
    "FOUW017": np.array([+5, -2]),
    "FOUW010": np.array([-1, +1]),

    "WSA101": np.array([+1, +1]),
    "WSA136": np.array([-1, -1]),

    "FOUW034": np.array([+1, 0]),
    "FOUW038": np.array([-1, 0]),

    "WSA124": np.array([-1, +1]),
    "WSA135": np.array([+1, -1])
}

# scale all nudges
for k in nudges:
    nudges[k] = nudges[k] * scale

snap_to_drainage(
    flow_dirs_filename="welsh_d8_fix_bounds.tif",
    sample_sites_filename="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Unique_SiteID_XY.csv",
    drainage_area_threshold=1000000,  # 1 km^2
    plot=True,
    save=True,
    nudges=nudges
    # nudges = {
    # "FOUW017": np.array([+1, -1]),   # east, south
    # "FOUW010": np.array([-1, +1]),    # west, north
    # "WSA101": np.array([+1, +1]),   # east, north
    # "WSA136": np.array([-1, -1]),    # west, south
    # "FOUW034": np.array([+1, 0]),   # east
    # "FOUW038": np.array([-1, 0]),    # west
    # "WSA124": np.array([-1, +1]),   # west, north
    # "WSA135": np.array([+1, -1])    # east, south
    # },
)

# Once this is done, we can load in the snapped sample sites and build the sample network again.
# Load sample network
sample_network, labels = get_sample_graph(
    flowdirs_filename="welsh_d8_fix_bounds.tif",
    sample_data_filename="Unique_SiteID_XY_snapped.csv",
)

plt.imshow(labels)
plt.axis('off')
plt.show()

plt.title("Correctly Connected Network with Snapped Samples")
plt.axis('off')
plot_network(sample_network)
plt.show()


# The network should now be connected properly, but should be checked to ensure that the snapping was successful and the samples have been snapped to the correct part of the network.

"""
Step 6: Creating CSV file with snapped coordinates and original coordinates
"""

# Step 1: Load original and snapped coordinates
original = pd.read_csv("/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Unique_SiteID_XY.csv")

# The snap_to_drainage() function should have created a new CSV.
snapped = pd.read_csv("Unique_SiteID_XY_snapped.csv")

print("Original columns:", original.columns)
print("Snapped columns:", snapped.columns)

# Step 2: Merge original + snapped data
# Rename snapped columns to distinguish them
snapped = snapped.rename(columns={
    "x": "snapped_x",
    "y": "snapped_y"
})

# Merge on the unique site identifier (adjust column name as needed)
merged = pd.merge(original, snapped[["Site ID", "snapped_x", "snapped_y"]], on="Site ID", how="left")

# Step 3: Save the merged dataset
merged.to_csv("/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Unique_SiteID_XY_snapped_original.csv", index=False)
print("✅ Saved combined CSV with original and snapped coordinates as 'Unique_SiteID_XY_snapped_original.csv'")

"""
Step 7: Add watercourse information to merged dataset with orginal and snapped x and y
"""

# Step 1: Load both datasets
snapped = pd.read_csv("/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Unique_SiteID_XY_snapped_original.csv")
phosphate = pd.read_csv("/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Filtered_Merged_Phosphate_Sites_XY.csv")

# Step 2: Reduce phosphate data to unique site-level info
phosphate_unique = phosphate[[
    "Site ID",
    "Watercourse (if known)",
    "Catchment (or \"tbc\")"
]].drop_duplicates(subset="Site ID")

# Step 3: Merge on Site ID
merged = pd.merge(
    snapped,
    phosphate_unique,
    on="Site ID",
    how="left"
)

# Step 4: Select desired columns and order
final = merged[[
    "Site ID",
    "x",
    "y",
    "snapped_x",
    "snapped_y",
    "Watercourse (if known)",
    "Catchment (or \"tbc\")"
]]

# Step 5: Save result
final.to_csv("/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Output/Snapped_Watercourse.csv", index=False)

print("✅ Saved merged site-level dataset as 'Snapped_Watercourse.csv'")
print(f"Rows in final file: {len(final)} (one per Site ID)")
