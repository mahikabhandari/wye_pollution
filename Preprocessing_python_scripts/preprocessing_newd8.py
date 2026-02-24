"""
This script preprocesses the input citizen science phosphate data from the Wye as well as a D8 flow direction raster so that they can be used for further analysis and inverse modeling.

Step 1: Merge the phosphate data with site information.
Step 2: Add x and y coordinates in BNG.
Step 3: Filter sites based on date and location.
Step 4: Check D8 raster is valid for use in inverse modelling software
Step 5: Prepare data for inverse modelling by snapping sites to D8 flow paths.
"""

#!/usr/bin/env python3
import os

import matplotlib.pyplot as plt
import pandas as pd
import rasterio
from pyproj import Transformer
from funmixer import (
    check_d8,
    get_sample_graph,
    plot_network,
    set_d8_boundaries_to_zero,
    snap_to_drainage,
)
###########################################################
# Step 1: Merge the phosphate data with site information.
###########################################################

# Load the two CSV files into dataframes
phosphate_df = pd.read_csv("Data_PP/Original/Phosphate_date.csv")
sites_df = pd.read_csv("Data_PP/Original/Sites.csv")

# Merge the dataframes on 'Site ID' column
merged_df = pd.merge(phosphate_df, sites_df, on='Site ID', how='left')

# Reorder columns to match the desired order
merged_df = merged_df[[
    'Sample ID', 'Date', 'Time', 'Phosphate (Hanna)', 'Site ID', 'Phosphate st', 'Hanna LR phos', 
    'Org', 'Name in EC', 'Latitude', 'Longitude', 'Watercourse (if known)', 'Catchment (or "tbc")'
]]

# Save the merged dataframe to a new CSV file
merged_df.to_csv("Data_PP/Merged_Phosphate_Sites.csv", index=False)

print("Merge complete. Saved as 'Merged_Phosphate_Sites.csv'")

###########################################################
# Step 2: Add x and y coordinates in BNG.
##########################################################


# Read the merged CSV
merged_df = pd.read_csv("Data_PP/Merged_Phosphate_Sites.csv")

#  Ensure column names are clean
merged_df.columns = merged_df.columns.str.strip()

# Define transformer for WGS84 -> British National Grid ---
transformer = Transformer.from_crs("EPSG:4326", "EPSG:27700", always_xy=True)

# Apply transformation to all rows
def latlon_to_xy(lon, lat):
    x, y = transformer.transform(lon, lat)
    return x, y

merged_df[['x', 'y']] = merged_df.apply(
    lambda row: latlon_to_xy(row['Longitude'], row['Latitude']),  
    axis=1,
    result_type='expand'
)

# Save the new CSV
merged_df.to_csv("Data_PP/Merged_Phosphate_Sites_XY.csv", index=False)

print("✅ Added x and y coordinates (British National Grid) and saved as 'Merged_Phosphate_Sites_XY.csv'")

############################################################
# Step 3: Filter sites

# Read the merged CSV
merged_df = pd.read_csv("Data_PP/Merged_Phosphate_Sites_XY.csv")

# Clean column names (optional but safe)
merged_df.columns = merged_df.columns.str.strip()

# Convert 'Date' column to datetime format
merged_df['Date'] = pd.to_datetime(merged_df['Date'], errors='coerce')

# Filter by date range (last two years: 2023-09-01 to 2025-09-01)
START_DATE = pd.Timestamp('2023-09-01')
END_DATE = pd.Timestamp('2025-09-01')
filtered_df = merged_df[(merged_df['Date'] >= START_DATE) & (merged_df['Date'] <= END_DATE)]

# Remove rows with missing or invalid Latitude/Longitude values ---
# Convert to numeric to handle cases where they're stored as strings
filtered_df['Latitude'] = pd.to_numeric(filtered_df['Latitude'], errors='coerce')
filtered_df['Longitude'] = pd.to_numeric(filtered_df['Longitude'], errors='coerce')

# Drop rows where Latitude or Longitude are missing (NaN)
filtered_df = filtered_df.dropna(subset=['Latitude', 'Longitude'])

# Remove rows where 'Hanna LR phos' is empty (NaN)
filtered_df = filtered_df.dropna(subset=['Hanna LR phos'])

SAMPLE_NUMBER_THRESHOLD = 50

# Filter Site IDs with at least 50 samples
site_counts = filtered_df['Site ID'].value_counts()
sites_with_50_samples = site_counts[site_counts >= SAMPLE_NUMBER_THRESHOLD].index
filtered_df = filtered_df[filtered_df['Site ID'].isin(sites_with_50_samples)]

# Save the filtered dataframe to a new CSV
output_path = "Data_PP/Filtered_Merged_Phosphate_Sites_XY.csv"
filtered_df.to_csv(output_path, index=False)

# Print summary info
unique_sites = filtered_df['Site ID'].nunique()

print(f"✅ Filtered citsci data saved as:\n{output_path}")
print(f"📊 Rows after filtering: {len(filtered_df)}")
print(f"📍 Unique sites after filtering: {unique_sites}")

# Load your original CSV
input_csv = "Data_PP/Filtered_Merged_Phosphate_Sites_XY.csv"
df = pd.read_csv(input_csv)

# Select only the desired columns
new_df = df[['Site ID', 'x', 'y']]

# Optional: ensure Sample ID is string
new_df['Site ID'] = new_df['Site ID'].astype(str)

# Drop duplicates based on Sample ID, keeping the first occurrence
unique_df = new_df.drop_duplicates(subset='Site ID', keep='first')

# Save to a new CSV
output_csv = "Data_PP/Unique_SiteID_XY.csv"
unique_df.to_csv(output_csv, index=False)

print(f"New CSV saved as {output_csv} with {len(unique_df)} unique rows.")

###########################################################
# Step 4: Preprocess the D8 raster
###########################################################

# Check if "Data_PP/Original/welsh_d8.nc" exists as a file and if not print a message saying it has to be downloaded manually
if not os.path.exists("Data_PP/Original/welsh_d8.nc"):
    print("File 'Data_PP/Original/welsh_d8.nc' not found. Please download it first (too big for repository!)")

# Check the validity of the file
check_d8("Data_PP/Original/welsh_d8.nc")

# Fix the boundary conditions using the set_d8_boundaries_to_zero function which sets all boundary cells to 0, writing the corrected raster to a new file.
set_d8_boundaries_to_zero("Data_PP/Original/welsh_d8.nc")
# set_d8_boundaries_to_zero("welsh_d8.nc")

# Now we can check the corrected raster and we confirm it is valid
fixed_raster_path = "welsh_d8_fix_bounds.tif"
check_d8(fixed_raster_path)


with rasterio.open(fixed_raster_path) as src:
    print("CRS:", src.crs)  # Coordinate Reference System
    print("Transform:", src.transform)

    # Pixel size in map units (usually meters if projected CRS)
    pixel_width = src.transform.a
    pixel_height = -src.transform.e  # usually negative due to raster orientation

    print(f"Pixel size: {pixel_width:.2f} x {pixel_height:.2f} (map units per pixel)")

###########################################################
# Step 5: Snap the sample sites to the drainage network
###########################################################

# In general, sample sites are not perfectly aligned with the drainage network, due to
# uncertain locations or the inherent simplification in representing flow using D8. This
# means that when genreating a "sample_network" using the get_sample_graphs function,
# the generated network may be incorrect. e.g., the sample sites may be connected to the wrong
# tributary or simply be disconnected from the network entirely. This can be fixed by snapping
# the sample sites to the nearest drainage network using the snap_to_drainage function.

# Load in real samples
samples = pd.read_csv("Data_PP/Unique_SiteID_XY.csv")
print(samples.columns)
sample_x, sample_y = samples["x"], samples["y"]

# Load sample network with unaligned samples
sample_network, labels = get_sample_graph(
    flowdirs_filename="welsh_d8_fix_bounds.tif",
    sample_data_filename="Data_PP/Unique_SiteID_XY.csv",
)

# plt.figure(figsize=(15, 10))  # Visualise network
plt.title("Disconnected Network Due to Misaligned Samples")
plot_network(sample_network)
plt.show()

# We can snap the noisy samples to the nearest drainage network using the snap_to_drainage function.
# This requires the flow direction raster and sample sites as input. The drainage_area_threshold parameter
# tells the function to only snap to drainage pixels with drainage area greater than the threshold. This is
# in the same units as the flow direction raster. The plot and save parameters control whether the snapped
# sample sites are plotted and saved to file respectively.

snap_to_drainage(
    flow_dirs_filename="welsh_d8_fix_bounds.tif",
    sample_sites_filename="Data_PP/Unique_SiteID_XY.csv",
    drainage_area_threshold=1000000,  # 1 km^2
    plot=True,
    save=True,
    #nudges={
    #"FOUW173": np.array([1000, 0])
    #}, 
    #nudges={"CG001": np.array([1000, -1000])},
    #nudges={
    #    "p1": np.array([1000, -1000]),
    #    "p6": np.array([500, 500]),
    #},
)

# Once this is done, we can load in the snapped sample sites and build the sample network again.
# Load sample network with snapped samples
sample_network, labels = get_sample_graph(
    flowdirs_filename="welsh_d8_fix_bounds.tif",
    sample_data_filename="Unique_SiteID_XY_snapped.csv",
)

plt.imshow(labels)
plt.axis('off')
plt.show()

# plt.figure(figsize=(15, 10))  # Visualise network
plt.title("Correctly Connected Network with Snapped Samples")
plot_network(sample_network)
plt.show()

# Load original and snapped coordinates ---
original = pd.read_csv("Data_PP/Unique_SiteID_XY.csv")

# The snap_to_drainage() function should have created a new CSV.
# Check what it's called — typically "Unique_SiteID_XY_snapped.csv"
snapped = pd.read_csv("Unique_SiteID_XY_snapped.csv")

# --- Step 3: Merge original + snapped data ---
# Rename snapped columns to distinguish them
snapped = snapped.rename(columns={
    "x": "snapped_x",
    "y": "snapped_y"
})

# Merge on the unique site identifier (adjust column name as needed)
merged = pd.merge(original, snapped[["Site ID", "snapped_x", "snapped_y"]], on="Site ID", how="left")

# Save the merged dataset ---
merged.to_csv("Data_PP/Unique_SiteID_XY_snapped_original.csv", index=False)
print("✅ Saved combined CSV with original and snapped coordinates as 'Unique_SiteID_XY_snapped_original.csv'")

# Add watercourse information to merged dataset with orginal and snapped x and y

# Load both datasets 
snapped = pd.read_csv("Data_PP/Unique_SiteID_XY_snapped_original.csv")
phosphate = pd.read_csv("Data_PP/Filtered_Merged_Phosphate_Sites_XY.csv")

# --- Reduce phosphate data to unique site-level info ---
# Keep only the relevant columns
phosphate_unique = phosphate[[
    "Site ID",
    "Watercourse (if known)",
    "Catchment (or \"tbc\")"
]].drop_duplicates(subset="Site ID")

# --- Merge on Site ID (site-level join) ---
merged = pd.merge(
    snapped,
    phosphate_unique,
    on="Site ID",
    how="left"
)

# --- Select desired columns and order ---
final = merged[[
    "Site ID",
    "x",
    "y",
    "snapped_x",
    "snapped_y",
    "Watercourse (if known)",
    "Catchment (or \"tbc\")"
]]

# --- Save result ---
# Check if the directory "Data_PP/Output" exists, if not, create it
if not os.path.exists("Data_PP/Output"):
    os.makedirs("Data_PP/Output")

final.to_csv("Data_PP/Output/Snapped_Watercourse.csv", index=False)

print("✅ Saved merged site-level dataset as 'Snapped_Watercourse.csv'")
print(f"Rows in final file: {len(final)} (one per Site ID)")
