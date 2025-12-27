import pandas as pd
import os
import matplotlib.pyplot as plt
import logging
import funmixer
logging.getLogger().addHandler(logging.StreamHandler())
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import rasterio
from matplotlib.colors import LogNorm

def find_duplicate_snapped_sites(csv_path="Snapped_Watercourse.csv"):
    """
    Reads a CSV containing site coordinate information and identifies all sites
    that share the same snapped_x and snapped_y coordinates.

    Required columns:
        - Site ID
        - x
        - y
        - snapped_x
        - snapped_y
        - Watercourse (if known)
        - Catchment (or "tbc")

    Returns a DataFrame containing only rows where two or more sites share
    the same snapped coordinates.
    """

    # Read CSV
    df = pd.read_csv(csv_path)

    # Ensure the needed columns exist
    required_columns = [
        "Site ID", "snapped_x", "snapped_y",
        "Watercourse (if known)", "Catchment (or \"tbc\")"
    ]
    for col in required_columns:
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")

    # Group by snapped coordinate pairs
    grouped = df.groupby(["snapped_x", "snapped_y"])

    # Filter groups with more than one site
    duplicates = grouped.filter(lambda g: len(g) > 1)

    return duplicates[
        ["Site ID", "snapped_x", "snapped_y",
         "Watercourse (if known)", "Catchment (or \"tbc\")"]
    ]

duplicates_df = find_duplicate_snapped_sites("/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Output/Snapped_Watercourse.csv")
print(duplicates_df)

def merge_csv(sites: str, samples: str, output_file: str):
    """
    Merge site CSV and sample CSV using 'Site ID'.
    
    Output columns:
    Sample ID, Date, Site ID, snapped_x, snapped_y, Phosphate (Hanna)
    """
    # Read the input CSVs
    sites = pd.read_csv(sites)
    samples = pd.read_csv(samples)

    # Merge on Site ID (inner join keeps only matching sites)
    merged = pd.merge(samples, sites, on='Site ID', how='inner')

    # Keep only required columns
    merged = merged[[
        "Sample ID",
        "Date",
        "Site ID",
        "snapped_x",
        "snapped_y",
        "Hanna LR phos"
    ]]

    # Save merged file
    merged.to_csv(output_file, index=False)
    print(f"Merged CSV saved to: {output_file}")
    print(f"\nNumber of rows in merged data: {len(merged)}")
    return merged

merge_csv(
    sites="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Output/Snapped_Watercourse.csv",
    samples="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Filtered_Merged_Phosphate_Sites_XY.csv",
    output_file="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Modelling/merged_output.csv"
)

"""
This script analyses a merged CSV of phosphate samples and site coordinates.
It computes:
1. Maximum and minimum Hanna LR phosphate values.
2. The Site ID with the most samples and the count.
3. The Site ID with the least samples and the count.
"""
# Load the CSV
df = pd.read_csv("/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Modelling/merged_output.csv")

# --- 1. Maximum and minimum Hanna LR phos ---
max_phos = df["Hanna LR phos"].max()
min_phos = df["Hanna LR phos"].min()

# --- 2. Site ID with most samples ---
site_counts = df["Site ID"].value_counts()
site_most_samples = site_counts.idxmax()
most_sample_count = site_counts.max()

# --- 3. Site ID with least samples ---
site_least_samples = site_counts.idxmin()
least_sample_count = site_counts.min()

# Print results
print("Maximum Hanna LR phos:", max_phos)
print("Minimum Hanna LR phos:", min_phos)

print("\nSite with MOST samples:", site_most_samples, 
      f"({most_sample_count} samples)")

print("Site with LEAST samples:", site_least_samples, 
      f"({least_sample_count} samples)")


def average_phosphate(
    merged_file: str,
    start_date: str,
    end_date: str,
    output_file: str
):
    """
    Computes the phosphate mean and standard deviation for each unique
    snapped_x / snapped_y location, automatically merging all Site IDs
    that share the same snapped coordinates.

    Output columns:
        Site IDs (merged list),
        snapped_x,
        snapped_y,
        Phosphate_mean,
        Phosphate_std
    """

    # Load data
    df = pd.read_csv(merged_file)

    # Convert date column
    df["Date"] = pd.to_datetime(df["Date"])

    # Filter by date range
    mask = (df["Date"] >= pd.to_datetime(start_date)) & (df["Date"] <= pd.to_datetime(end_date))
    df = df.loc[mask]

    # --- Merge sites by snapped_x / snapped_y ---
    # Combine all Site IDs per snapped coordinate
    site_groups = (
        df.groupby(["snapped_x", "snapped_y"])["Site ID"]
        .apply(lambda x: ", ".join(sorted(set(map(str, x)))))
        .reset_index()
        .rename(columns={"Site ID": "Site IDs"})
    )

    # --- Compute phosphate stats for each snapped coordinate pair ---
    stats = (
        df.groupby(["snapped_x", "snapped_y"])
        .agg(
            Phosphate_mean=("Hanna LR phos", "mean"),
            Phosphate_std=("Hanna LR phos", "std")
        )
        .reset_index()
    )

    # Merge site ID information with stats
      # Merge site ID information with stats
    result = pd.merge(site_groups, stats, on=["snapped_x", "snapped_y"])

    # Reorder columns so funmixer works
    result = result[
        ["Site IDs", "snapped_x", "snapped_y",
         "Phosphate_mean", "Phosphate_std"]
    ]

    # Save CSV
    result.to_csv(output_file, index=False)
    print(f"Averaged phosphate CSV saved to: {output_file}")
    print(f"Number of snapped coordinate groups: {len(result)}")

    return result

average_phosphate(
    merged_file="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Modelling/merged_output.csv",
    start_date="2023-10-19",
    end_date="2025-10-19",
    output_file="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Modelling/phosphate_averages.csv"
)

"""
This script is a minimum working example for how to unmix the downstream observations for a specific element.
It treats each sub-basin defined by a sample as a discrete variable, and solves for the upstream concentrations of the element in each sub-basin.

It loads a sample network graph and observations from data files, visualizes the network, and builds the optimization problem.
Then, it performs a sweep of different regularization strengths for the problem and visualizes the results.
Next, it solves the problem using a specified solver and regularization strength, obtaining predicted downstream and upstream concentrations.
The script also calculates unique upstream areas for each basin in the network and generates an upstream concentration map based on the predictions.
Finally, it visualizes the predicted downstream concentrations and the upstream concentration map for the specified element (default: Mg).
"""

def run_funmixer(unmix_csv: str, flowdir_file: str, element_col: str):
    """
    Run funmixer network unmixing on the provided CSV and flow direction raster.
    """
    obs_data = pd.read_csv(unmix_csv)

    # Keep only required columns for funmixer
    required_cols = ["Site IDs", element_col]
    obs_data = obs_data[required_cols]

    # Load network
    sample_network, labels = funmixer.get_sample_graph(
        flowdirs_filename=flowdir_file,
        sample_data_filename=unmix_csv
    )

    plt.figure(figsize=(15, 10))
    funmixer.plot_network(sample_network)
    plt.show()
    print("Building funmixer problem...")
    problem = funmixer.SampleNetworkUnmixer(sample_network=sample_network)
    element_data = funmixer.get_element_obs(element_col, obs_data)

    print(element_data.keys())

    nodes = set(list(sample_network.nodes))
    obsnodes = set(list(element_data.keys()))
    missing = nodes - obsnodes
    # Check if set sare the same
    same = nodes == obsnodes

    # Print length of each set
    print(f"Number of nodes in network: {len(nodes)}")
    print(f"Number of observation nodes: {len(obsnodes)}")

    print(f"Sets are the same: {same}")
    print(f"Missing observation nodes: {missing}")

    funmixer.plot_sweep_of_regularizer_strength(problem, element_data, -5, -1, 11)
    regularizer_strength = 10 ** (-3.4)
    print(f"Chosen regularization strength: {regularizer_strength}")

    print("Solving problem...")
    solution = problem.solve(element_data, solver="clarabel", regularization_strength=regularizer_strength)

    area_dict = funmixer.get_unique_upstream_areas(sample_network, labels)
    upstream_map = funmixer.get_upstream_concentration_map(area_dict, solution.upstream_preds)

    # Visualise downstream
    funmixer.visualise_downstream(pred_dict=solution.downstream_preds, obs_dict=element_data, element=element_col)
    plt.show()

    # Visualise upstream
    plt.imshow(upstream_map)
    cb = plt.colorbar()
    cb.set_label(f"{element_col} concentration")
    plt.title("Upstream Concentration Map")
    plt.show()

    # Load raster metadata from the flow direction file
    with rasterio.open(flowdir_file) as src:
        transform = src.transform
        crs = src.crs
        bounds = src.bounds
        extent = [bounds.left, bounds.right, bounds.bottom, bounds.top]

    fig = plt.figure(figsize=(10, 10))
    ax = plt.axes(projection=ccrs.OSGB())

    ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax.add_feature(cfeature.RIVERS, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linestyle="--", alpha=0.5)
    ax.gridlines(draw_labels=True)

    im = ax.imshow(
        upstream_map,
        origin="upper",
        extent=extent,
        transform=ccrs.OSGB(),
        cmap="viridis",
        norm=LogNorm(),
        interpolation="nearest"
    )

    ax.set_extent([
        bounds.left + (bounds.right - bounds.left) * 0.2,
        bounds.right - (bounds.right - bounds.left) * 0.2,
        bounds.bottom + (bounds.top - bounds.bottom) * 0.2,
        bounds.top - (bounds.top - bounds.bottom) * 0.2,
    ], crs=ccrs.OSGB())

    cb = plt.colorbar(im, ax=ax, shrink=0.7)
    cb.set_label(f"{element_col} concentration (log scale)")

    plt.title(f"Upstream Concentration Map for {element_col}")
    plt.show()

run_funmixer(
    unmix_csv="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Modelling/phosphate_averages.csv",
    flowdir_file="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Original/welsh_d8.nc",
    element_col="Phosphate_mean"
)
