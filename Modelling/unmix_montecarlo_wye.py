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
import matplotlib.colors as mcolors
from matplotlib.ticker import FormatStrFormatter
import numpy as np

""" Sorting data before applying funmixer unmixing """

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
Monte Carlo version of funmixer unmixing to propagate observational uncertainty WITHOUT RIVER NETWORK

This script:
- Loads downstream observations and a flow direction raster
- Builds a sample network
- Solves the unmixing problem using Monte Carlo resampling
- Computes mean and uncertainty of downstream and upstream predictions
- Saves misfit and uncertainty to CSV
- Visualises upstream concentration and uncertainty maps
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import rasterio
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import LogNorm
from matplotlib.ticker import FormatStrFormatter
import funmixer


def run_funmixer_montecarlo(
    unmix_csv: str,
    flowdir_file: str,
    element_col: str,
    relative_error: float = 10.0,   # % uncertainty on observations
    num_repeats: int = 50,
    regularization_strength: float = 10 ** (-3.4),
):
    """
    Run funmixer network unmixing with Monte Carlo uncertainty propagation.
    """

    # ------------------------------------------------------------------
    # Load observations
    # ------------------------------------------------------------------
    obs_data = pd.read_csv(unmix_csv)
    obs_data = obs_data[["Site IDs", element_col]]

    # ------------------------------------------------------------------
    # Load network
    # ------------------------------------------------------------------
    sample_network, labels = funmixer.get_sample_graph(
        flowdirs_filename=flowdir_file,
        sample_data_filename=unmix_csv,
    )

    plt.figure(figsize=(15, 10))
    funmixer.plot_network(sample_network)
    plt.title("Sample Network")
    plt.show()

    print("Building funmixer problem...")
    problem = funmixer.SampleNetworkUnmixer(sample_network=sample_network)

    element_data = funmixer.get_element_obs(element_col, obs_data)
    # --------------------------------------------------
    # Regularisation sweep (diagnostic)
    # --------------------------------------------------
    plt.figure(figsize=(6, 4))
    funmixer.plot_sweep_of_regularizer_strength(
        problem,
        element_data,
        -5,   # 10^-5
        -1,   # 10^-1
        11,   # number of values
    )
    plt.tight_layout()
    plt.savefig(
        "/Users/mahikabhandari/Desktop/Earth Science Year 4/"
        "MSci Project/Analysis/Modelling/Figures/regularisation_sweep_nrn.png",
        dpi=300,
    )
    plt.show()

    # Choose regularisation strength
    regularization_strength = 10 ** (-3.4)
    print(f"Chosen regularization strength: {regularization_strength}")

    # --------------------------------------------------
    # Deterministic solve (diagnostic only)
    # --------------------------------------------------
    solution = problem.solve(
        element_data,
        solver="clarabel",
        regularization_strength=regularization_strength,
    )

    # Observed vs predicted downstream
    plt.figure(figsize=(5, 5))
    funmixer.visualise_downstream(
        pred_dict=solution.downstream_preds,
        obs_dict=element_data,
        element=element_col,
    )
    plt.tight_layout()
    plt.savefig(
        "/Users/mahikabhandari/Desktop/Earth Science Year 4/"
        "MSci Project/Analysis/Modelling/Figures/predicted_vs_observed_nrn.png",
        dpi=300,
    )
    plt.show()

    # ------------------------------------------------------------------
    # Monte Carlo solve
    # ------------------------------------------------------------------
    print("Running Monte Carlo unmixing...")
    downstream_mc, upstream_mc = problem.solve_montecarlo(
        element_data,
        relative_error=relative_error,
        num_repeats=num_repeats,
        regularization_strength=regularization_strength,
    )

    # ------------------------------------------------------------------
    # Compute means and uncertainties
    # ------------------------------------------------------------------
    downstream_mean = {k: np.mean(v) for k, v in downstream_mc.items()}
    downstream_std = {k: np.std(v) for k, v in downstream_mc.items()}

    upstream_mean = {k: np.mean(v) for k, v in upstream_mc.items()}
    upstream_std = {k: np.std(v) for k, v in upstream_mc.items()}

    # ------------------------------------------------------------------
    # Misfit + uncertainty table
    # ------------------------------------------------------------------
    misfit_df = pd.DataFrame({
        "Site_IDs": list(element_data.keys()),
        "Observed": list(element_data.values()),
        "Predicted_mean": [downstream_mean[k] for k in element_data.keys()],
        "Predicted_std": [downstream_std[k] for k in element_data.keys()],
    })

    misfit_df["Log_Misfit"] = np.log(
        misfit_df["Observed"] / misfit_df["Predicted_mean"]
    )

    misfit_file = (
        "/Users/mahikabhandari/Desktop/Earth Science Year 4/"
        "MSci Project/Analysis/Modelling/phosphate_montecarlo_misfit.csv"
    )
    misfit_df.to_csv(misfit_file, index=False)
    print(f"Monte Carlo misfit CSV saved to:\n{misfit_file}")

    # ------------------------------------------------------------------
    # Build upstream maps
    # ------------------------------------------------------------------
    area_dict = funmixer.get_unique_upstream_areas(sample_network, labels)

    upstream_mean_map = funmixer.get_upstream_concentration_map(
        area_dict, upstream_mean
    )

    upstream_std_map = funmixer.get_upstream_concentration_map(
        area_dict, upstream_std
    )

    # ------------------------------------------------------------------
    # Load raster metadata
    # ------------------------------------------------------------------
    with rasterio.open(flowdir_file) as src:
        bounds = src.bounds
        extent = [bounds.left, bounds.right, bounds.bottom, bounds.top]

    # Wye catchment bounds (OSGB, meters)
    wye_bounds = [300000, 350000, 220000, 270000]  # [west, east, south, north]

    # Add a margin to zoom out
    margin = 30000  # 30 km
    wye_bounds_zoomed_out = [
        wye_bounds[0] - margin,
        wye_bounds[1] + margin,
        wye_bounds[2] - margin,
        wye_bounds[3] + margin,
    ]

    # ------------------------------------------------------------------
    # Plot upstream mean
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={"projection": ccrs.OSGB()})
    ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax.add_feature(cfeature.RIVERS, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linestyle="--", alpha=0.5)

    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={"projection": ccrs.OSGB()})
    ax.set_extent(wye_bounds_zoomed_out, crs=ccrs.OSGB())
        # -------------------------------
    # Axis ticks (OSGB coordinates)
    # -------------------------------
    xmin, xmax, ymin, ymax = ax.get_extent(crs=ccrs.OSGB())

    x_ticks = np.arange(
        np.floor(xmin / 20_000) * 20_000,
        xmax + 1,
        20_000,
    )
    y_ticks = np.arange(
        np.floor(ymin / 20_000) * 20_000,
        ymax + 1,
        20_000,
    )

    ax.set_xticks(x_ticks, crs=ccrs.OSGB())
    ax.set_yticks(y_ticks, crs=ccrs.OSGB())

    ax.xaxis.set_major_formatter(FormatStrFormatter("%d"))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%d"))

    ax.set_xlabel("Easting (m, OSGB)")
    ax.set_ylabel("Northing (m, OSGB)")


    im = ax.imshow(
        upstream_mean_map,
        origin="upper",
        extent=extent,
        transform=ccrs.OSGB(),
        cmap="viridis",
        norm=LogNorm(),
    )

    cb = plt.colorbar(im, ax=ax, shrink=0.7)
    cb.set_label(f"{element_col} mean concentration (log scale)")
    plt.title(f"Upstream Mean Concentration ({element_col})")
    plt.show()

    # ------------------------------------------------------------------
    # Plot upstream uncertainty
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={"projection": ccrs.OSGB()})
    ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax.add_feature(cfeature.RIVERS, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linestyle="--", alpha=0.5)

    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={"projection": ccrs.OSGB()})
    ax.set_extent(wye_bounds_zoomed_out, crs=ccrs.OSGB())
    # -------------------------------
    # Axis ticks (OSGB coordinates)
    # -------------------------------
    xmin, xmax, ymin, ymax = ax.get_extent(crs=ccrs.OSGB())

    x_ticks = np.arange(
        np.floor(xmin / 20_000) * 20_000,
        xmax + 1,
        20_000,
    )
    y_ticks = np.arange(
        np.floor(ymin / 20_000) * 20_000,
        ymax + 1,
        20_000,
    )

    ax.set_xticks(x_ticks, crs=ccrs.OSGB())
    ax.set_yticks(y_ticks, crs=ccrs.OSGB())

    ax.xaxis.set_major_formatter(FormatStrFormatter("%d"))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%d"))

    ax.set_xlabel("Easting (m, OSGB)")
    ax.set_ylabel("Northing (m, OSGB)")

    im = ax.imshow(
        upstream_std_map,
        origin="upper",
        extent=extent,
        transform=ccrs.OSGB(),
        cmap="magma",
        norm=LogNorm(),
    )

    cb = plt.colorbar(im, ax=ax, shrink=0.7)
    cb.set_label(f"{element_col} uncertainty (log scale)")
    plt.title(f"Upstream Uncertainty ({element_col})")
    plt.show()

    return downstream_mean, downstream_std, upstream_mean, upstream_std


# ----------------------------------------------------------------------
# Run
# ----------------------------------------------------------------------

run_funmixer_montecarlo(
    unmix_csv="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Modelling/phosphate_averages.csv",
    flowdir_file="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Original/welsh_d8.nc",
    element_col="Phosphate_mean",
    relative_error=10.0,
    num_repeats=50,
) 


"""
Monte Carlo version of funmixer unmixing to propagate observational uncertainty WITH RIVER NETWORK

This script:
- Loads downstream observations and a flow direction raster
- Builds a sample network
- Solves the unmixing problem using Monte Carlo resampling
- Computes mean and uncertainty of downstream and upstream predictions
- Saves misfit and uncertainty to CSV
- Visualises upstream concentration and uncertainty maps
  with the river network included in the figure
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import rasterio
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import LogNorm
from osgeo import gdal

import funmixer
import funmixer.flow_acc_cfuncs as cf

# ---------------------------------------------------------------------
# D8 accumulator (unchanged logic, simplified placement)
# ---------------------------------------------------------------------
class D8Accumulator:
    def __init__(self, filename: str):
        self.arr, self.ds = self._read(filename)
        self.arr = self.arr.astype(int)
        self.receivers = cf.d8_to_receivers(self.arr)
        base = np.where(self.receivers == np.arange(len(self.receivers)))[0]
        self.order = cf.build_ordered_list_iterative(self.receivers, base)

    def _read(self, filename):
        ds = gdal.Open(filename)
        if ds is None:
            raise IOError(f"Could not open {filename}")
        return ds.GetRasterBand(1).ReadAsArray(), ds

    def accumulate(self):
        return cf.accumulate_flow(
            self.receivers,
            self.order,
            weights=np.ones(len(self.receivers)),
        ).reshape(self.arr.shape)

    def indices_to_coords(self, rows, cols):
        gt = self.ds.GetGeoTransform()
        x = gt[0] + cols * gt[1] + gt[1] / 2
        y = gt[3] + rows * gt[5] + gt[5] / 2
        return x, y

# ---------------------------------------------------------------------
# River network builder
# ---------------------------------------------------------------------
def build_river_network(flowdir_file, drainage_area_threshold=5e6):
    accum = D8Accumulator(flowdir_file)

    gt = accum.ds.GetGeoTransform()
    cell_area = abs(gt[1] * gt[5])

    area = accum.accumulate() * cell_area
    rows, cols = np.where(area > drainage_area_threshold)

    x, y = accum.indices_to_coords(rows, cols)
    loga = np.log10(area[rows, cols])
    size = 6 * (loga - loga.min() + 0.05) / (loga.max() - loga.min())

    return x, y, size

# ---------------------------------------------------------------------
# Plot upstream map + river network
# ---------------------------------------------------------------------
def plot_upstream_with_rivers(
    upstream_map,
    extent,
    river_x,
    river_y,
    river_size,
    title,
    cbar_label,
    cmap,
    output_path,
):
    fig, ax = plt.subplots(
        figsize=(10, 10),
        subplot_kw={"projection": ccrs.OSGB()},
    )

    ax.add_feature(cfeature.COASTLINE, linewidth=0.8, zorder=0)
    ax.add_feature(cfeature.BORDERS, linestyle=":", linewidth=0.8, zorder=0)

    ax.scatter(
        river_x,
        river_y,
        s=river_size,
        c="steelblue",
        alpha=0.5,
        transform=ccrs.OSGB(),
        zorder=1,
    )

    im = ax.imshow(
        upstream_map,
        origin="upper",
        extent=extent,
        transform=ccrs.OSGB(),
        cmap=cmap,
        norm=LogNorm(),
        zorder=2,
    )

    wye = [300000, 350000, 220000, 270000]
    margin = 30000
    ax.set_extent(
        [wye[0]-margin, wye[1]+margin, wye[2]-margin, wye[3]+margin],
        crs=ccrs.OSGB(),
    )
    # -------------------------------
    # Axis ticks (OSGB coordinates)
    # -------------------------------
    xmin, xmax, ymin, ymax = ax.get_extent(crs=ccrs.OSGB())

    x_ticks = np.arange(
        np.floor(xmin / 20_000) * 20_000,
        xmax + 1,
        20_000,
    )
    y_ticks = np.arange(
        np.floor(ymin / 20_000) * 20_000,
        ymax + 1,
        20_000,
    )

    ax.set_xticks(x_ticks, crs=ccrs.OSGB())
    ax.set_yticks(y_ticks, crs=ccrs.OSGB())

    ax.xaxis.set_major_formatter(FormatStrFormatter("%d"))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%d"))

    ax.set_xlabel("Easting (m, OSGB)")
    ax.set_ylabel("Northing (m, OSGB)")

    cbar = plt.colorbar(im, ax=ax, shrink=0.7)
    cbar.set_label(cbar_label)

    ax.set_title(title)
    ax.set_aspect("equal")

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

# ---------------------------------------------------------------------
# Main Monte Carlo workflow
# ---------------------------------------------------------------------
    
def run_funmixer_montecarlo(
    unmix_csv: str,
    flowdir_file: str,
    element_col: str,
    relative_error: float = 10.0,   # % uncertainty on observations
    num_repeats: int = 50,
    regularization_strength: float = 10 ** (-3.4),
    drainage_area_threshold=5e6,
):
    # -------------------------------
    # Load observations + network
    # -------------------------------
    obs = pd.read_csv(unmix_csv)[["Site IDs", element_col]]

    sample_network, labels = funmixer.get_sample_graph(
        flowdirs_filename=flowdir_file,
        sample_data_filename=unmix_csv,
    )

    problem = funmixer.SampleNetworkUnmixer(sample_network)
    element_data = funmixer.get_element_obs(element_col, obs)
    # --------------------------------------------------

    # Regularisation sweep (diagnostic)
    # --------------------------------------------------
    plt.figure(figsize=(6, 4))
    funmixer.plot_sweep_of_regularizer_strength(
        problem,
        element_data,
        -5,   # 10^-5
        -1,   # 10^-1
        11,   # number of values
    )
    plt.tight_layout()
    plt.savefig(
        "/Users/mahikabhandari/Desktop/Earth Science Year 4/"
        "MSci Project/Analysis/Modelling/Figures/regularisation_sweep.png",
        dpi=300,
    )
    plt.show()

    # Choose regularisation strength
    regularization_strength = 10 ** (-3.4)
    print(f"Chosen regularization strength: {regularization_strength}")

    # --------------------------------------------------
    # Deterministic solve (diagnostic only)
    # --------------------------------------------------
    solution = problem.solve(
        element_data,
        solver="clarabel",
        regularization_strength=regularization_strength,
    )

    # Observed vs predicted downstream
    plt.figure(figsize=(5, 5))
    funmixer.visualise_downstream(
        pred_dict=solution.downstream_preds,
        obs_dict=element_data,
        element=element_col,
    )
    plt.tight_layout()
    plt.savefig(
        "/Users/mahikabhandari/Desktop/Earth Science Year 4/"
        "MSci Project/Analysis/Modelling/Figures/predicted_vs_observed.png",
        dpi=300,
    )
    plt.show()


    # -------------------------------
    # Monte Carlo solve
    # -------------------------------

    downstream_mc, upstream_mc = problem.solve_montecarlo(
    element_data,
    relative_error=relative_error,
    num_repeats=num_repeats,
    regularization_strength=regularization_strength,
    )

    downstream_mean = {k: np.mean(v) for k, v in downstream_mc.items()}
    downstream_std  = {k: np.std(v)  for k, v in downstream_mc.items()}

    upstream_mean = {k: np.mean(v) for k, v in upstream_mc.items()}
    upstream_std  = {k: np.std(v)  for k, v in upstream_mc.items()}

    # -------------------------------
    # Save misfit CSV
    # -------------------------------
    misfit_df = pd.DataFrame({
        "Site_IDs": list(element_data.keys()),
        "Observed": list(element_data.values()),
        "Predicted_mean": [downstream_mean[k] for k in element_data],
        "Predicted_std":  [downstream_std[k]  for k in element_data],
    })

    misfit_df["Log_Misfit"] = np.log(
        misfit_df["Observed"] / misfit_df["Predicted_mean"]
    )

    out_csv = "phosphate_montecarlo_misfit.csv"
    misfit_df.to_csv(out_csv, index=False)
    print(f"Saved misfit CSV → {out_csv}")

    # -------------------------------
    # Upstream maps
    # -------------------------------
    area_dict = funmixer.get_unique_upstream_areas(sample_network, labels)

    upstream_mean_map = funmixer.get_upstream_concentration_map(
        area_dict, upstream_mean
    )
    upstream_std_map = funmixer.get_upstream_concentration_map(
        area_dict, upstream_std
    )

    with rasterio.open(flowdir_file) as src:
        extent = [src.bounds.left, src.bounds.right,
                  src.bounds.bottom, src.bounds.top]

    # -------------------------------
    # River network (once)
    # -------------------------------
    rx, ry, rs = build_river_network(
        flowdir_file,
        drainage_area_threshold,
    )

    # -------------------------------
    # Plots
    # -------------------------------
    plot_upstream_with_rivers(
        upstream_mean_map,
        extent,
        rx, ry, rs,
        title="Upstream Mean Phosphate Concentration",
        cbar_label="Mean concentration (log scale)",
        cmap="viridis",
        output_path="Modelling/Monte Carlo Figures/upstream_mean_montecarlo.png",
    )

    plot_upstream_with_rivers(
        upstream_std_map,
        extent,
        rx, ry, rs,
        title="Upstream Phosphate Uncertainty",
        cbar_label="Uncertainty (log scale)",
        cmap="magma",
        output_path="Modelling/Monte Carlo Figures/upstream_uncertainty_montecarlo.png",
    )

# ---------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------
run_funmixer_montecarlo(
    unmix_csv="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Modelling/phosphate_averages.csv",
    flowdir_file="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Original/welsh_d8.nc",
    element_col="Phosphate_mean",
)
