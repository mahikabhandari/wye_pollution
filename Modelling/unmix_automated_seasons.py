"""
Monte Carlo version of funmixer unmixing to propagate observational uncertainty WITH RIVER NETWORK.

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

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import funmixer
import funmixer.flow_acc_cfuncs as cf
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from matplotlib.colors import LogNorm
from matplotlib.ticker import FormatStrFormatter
from osgeo import gdal
from scipy.stats import gmean as geo_mean

PHOSPHATE_DECAY_K = 1.45e-6  # m^-1 - Universal phosphate decay constant


# ---------------------------------------------------------------------
# D8 accumulator (unchanged logic, simplified placement)
# ---------------------------------------------------------------------
class D8Accumulator:
    """D8 flow direction accumulator (borrowed and adapted from funmixer)."""

    def __init__(self, filename: str):
        """Initialize the D8Accumulator with a flow direction raster file."""
        self.arr, self.ds = self._read(filename)
        self.arr = self.arr.astype(int)
        self.receivers = cf.d8_to_receivers(self.arr)
        base = np.where(self.receivers == np.arange(len(self.receivers)))[0]
        self.order = cf.build_ordered_list_iterative(self.receivers, base)

    def _read(self, filename):
        ds = gdal.Open(filename)
        if ds is None:
            raise OSError(f"Could not open {filename}")
        return ds.GetRasterBand(1).ReadAsArray(), ds

    def accumulate(self):
        """Accumulate flow using the D8 flow direction."""
        return cf.accumulate_flow(
            self.receivers,
            self.order,
            weights=np.ones(len(self.receivers)),
        ).reshape(self.arr.shape)

    def indices_to_coords(self, rows, cols):
        """Convert raster indices to geographic coordinates."""
        gt = self.ds.GetGeoTransform()
        x = gt[0] + cols * gt[1] + gt[1] / 2
        y = gt[3] + rows * gt[5] + gt[5] / 2
        return x, y


# ---------------------------------------------------------------------
# River network builder
# ---------------------------------------------------------------------
def build_river_network(flowdir_file, drainage_area_threshold=5e6):
    """Build a river network (xs, ys and areas) from a flow direction raster."""
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
    cmap,
    vmin,
    vmax,
    ax=None,
):
    """Plot upstream concentration map with river network."""
    if ax is None:
        fig, ax = plt.subplots(
            figsize=(6, 6),
            subplot_kw={"projection": ccrs.OSGB()},
        )

    ax.add_feature(cfeature.COASTLINE, linewidth=0.5, zorder=0)
    ax.add_feature(cfeature.BORDERS, linestyle=":", linewidth=0.5, zorder=0)

    ax.scatter(
        river_x,
        river_y,
        s=river_size,
        c="steelblue",
        alpha=0.4,
        transform=ccrs.OSGB(),
        zorder=4,
    )

    im = ax.imshow(
        upstream_map,
        origin="upper",
        extent=extent,
        transform=ccrs.OSGB(),
        cmap=cmap,
        norm=LogNorm(vmin=vmin, vmax=vmax),
        zorder=2,
    )

    ax.set_title(title, fontsize=10)
    ax.set_aspect("equal")

    wye = [270000, 370000, 190000, 290000]
    margin = 30000
    ax.set_extent(
        [wye[0] - margin, wye[1] + margin, wye[2] - margin, wye[3] + margin],
        crs=ccrs.OSGB(),
    )
    # -------------------------------
    # Axis ticks (OSGB coordinates)
    # -------------------------------
    xmin, xmax, ymin, ymax = ax.get_extent(crs=ccrs.OSGB())

    x_ticks = np.arange(
        np.floor(xmin / 40_000) * 40_000,
        xmax + 1,
        40_000,
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

    return im


def format_season_title(season_key: str) -> str:
    """Format the season title for display."""
    parts = season_key.split("_")

    if len(parts) == 2:
        return f"{parts[0]} {parts[1]}"
    elif len(parts) == 3:
        return f"{parts[0]} {parts[1]}–{parts[2]}"
    else:
        return season_key.replace("_", " ")


# ---------------------------------------------------------------------
# Main Monte Carlo workflow
# ---------------------------------------------------------------------


def run_funmixer_montecarlo(
    unmix_csv: str,
    runoff_csv: str,
    flowdir_file: str,
    element_col: str,
    season: str,
    output_base: str,
    relative_error: float = 10.0,
    num_repeats: int = 50,
    regularization_strength: float = 10 ** (-3.4),
    drainage_area_threshold=5e6,
):
    """Run Monte Carlo version of funmixer unmixing to propagate observational uncertainty WITH RIVER NETWORK."""
    title_string = format_season_title(season)
    # -------------------------------
    # Load observations + network
    # -------------------------------
    obs = pd.read_csv(unmix_csv)[["Site_IDs", element_col]]

    sample_network, labels = funmixer.get_sample_graph(
        flowdirs_filename=flowdir_file,
        sample_data_filename=unmix_csv,
    )
    # print("Graph nodes:", list(sample_network.nodes))
    element_data = funmixer.get_element_obs(element_col, obs)

    # -------------------------------
    # Define rate constants
    # -------------------------------

    rate_constants = {site_id: PHOSPHATE_DECAY_K for site_id in element_data.keys()}

    print(
        "Decay constant range:",
        min(rate_constants.values()),
        max(rate_constants.values()),
    )

    print("Element data sites:", len(element_data))
    print("Network sites:", len(sample_network.nodes))
    print(
        "Missing in rate_constants:",
        set(sample_network.nodes) - set(rate_constants.keys()),
    )

    def load_runoff_export_rates(
        runoff_csv: str,
        site_col: str = "Site_IDs",
        runoff_col: str = "Mean_Runoff",
        normalise: bool = True,
        min_value: float = 1e-6,
    ) -> dict[str, float]:
        """Load runoff CSV (already merged by snapped coordinates) and return export rates for funmixer."""
        df = pd.read_csv(runoff_csv)

        export_rates = {
            str(row[site_col]): max(float(row[runoff_col]), min_value)
            for _, row in df.iterrows()
        }

        if normalise:
            gm = geo_mean(list(export_rates.values()))
            export_rates = {k: v / gm for k, v in export_rates.items()}

        return export_rates

    runoff_rates = load_runoff_export_rates(runoff_csv=runoff_csv)

    runoff_rates = {k: v for k, v in runoff_rates.items() if k in sample_network.nodes}

    missing = set(element_data.keys()) - set(runoff_rates.keys())
    extra = set(runoff_rates.keys()) - set(element_data.keys())

    print("Missing runoff for sites:", missing)
    print("Extra runoff sites:", extra)

    problem = funmixer.SampleNetworkUnmixer(
        sample_network,
        rate_constants=rate_constants,
    )

    # --------------------------------------------------

    # Choose regularisation strength
    regularization_strength = 10 ** (-3.4)
    print(f"Chosen regularization strength: {regularization_strength}")

    # --------------------------------------------------
    # Deterministic solve
    # --------------------------------------------------
    solution = problem.solve(
        element_data,
        solver="clarabel",
        regularization_strength=regularization_strength,
        export_rates=runoff_rates,
    )

    # Observed vs predicted downstream
    plt.figure(figsize=(7, 5))
    funmixer.visualise_downstream(
        pred_dict=solution.downstream_preds,
        obs_dict=element_data,
        element=element_col,
    )
    # Labels and title
    plt.xlabel("Observed Phosphate (mg/L)", fontsize=12)
    plt.ylabel("Predicted Phosphate (mg/L)", fontsize=12)
    plt.title("Predicted vs Observed Phosphate Concentration", fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(output_base, f"predicted_vs_observed_{season}.png"))

    # -------------------------------
    # Monte Carlo solve
    # -------------------------------

    downstream_mc, upstream_mc = problem.solve_montecarlo(
        element_data,
        relative_error=relative_error,
        num_repeats=num_repeats,
        export_rates=runoff_rates,
        regularization_strength=regularization_strength,
    )

    downstream_mean = {k: np.mean(v) for k, v in downstream_mc.items()}
    downstream_std = {k: np.std(v) for k, v in downstream_mc.items()}

    upstream_mean = {k: np.mean(v) for k, v in upstream_mc.items()}
    upstream_std = {k: np.std(v) for k, v in upstream_mc.items()}

    # --------------------------------
    # Compute Monte Carlo flux per site
    # --------------------------------
    conversion_factor = 1e-6 * 4

    upstream_flux_mean = {
        site: runoff_rates[site] * upstream_mean[site] * conversion_factor
        for site in upstream_mean
    }

    upstream_flux_std = {
        site: runoff_rates[site] * upstream_std[site] * conversion_factor
        for site in upstream_std
    }

    # --------------------------------
    # Save flux mean + uncertainty CSV
    # --------------------------------
    flux_df = pd.DataFrame(
        {
            "Site_IDs": list(upstream_flux_mean.keys()),
            "Phosphate_Flux_Mean_kg_m2_day": list(upstream_flux_mean.values()),
            "Phosphate_Flux_Std_kg_m2_day": [
                upstream_flux_std[site] for site in upstream_flux_mean
            ],
        }
    )

    flux_df["Season"] = season

    flux_csv = os.path.join(output_base, f"phosphate_flux_montecarlo_{season}.csv")

    flux_df.to_csv(flux_csv, index=False)
    print(f"Saved Monte Carlo flux CSV → {flux_csv}")

    print(
        "Downstream mean (with decay) min/max:",
        min(downstream_mean.values()),
        max(downstream_mean.values()),
    )

    # -------------------------------
    # Save misfit CSV
    # -------------------------------
    misfit_df = pd.DataFrame(
        {
            "Site_IDs": list(element_data.keys()),
            "Observed": list(element_data.values()),
            "Predicted_mean": [downstream_mean[k] for k in element_data],
            "Predicted_std": [downstream_std[k] for k in element_data],
        }
    )

    misfit_df["Log_Misfit"] = np.log(
        misfit_df["Observed"] / misfit_df["Predicted_mean"]
    )

    misfit_csv = os.path.join(output_base, f"phosphate_montecarlo_misfit_{season}.csv")

    misfit_df.to_csv(misfit_csv, index=False)
    print(f"Saved misfit CSV → {misfit_csv}")

    # -------------------------------
    # Upstream maps
    # -------------------------------
    area_dict = funmixer.get_unique_upstream_areas(sample_network, labels)

    upstream_mean_map = funmixer.get_upstream_concentration_map(
        area_dict, upstream_mean
    )
    upstream_std_map = funmixer.get_upstream_concentration_map(area_dict, upstream_std)

    with rasterio.open(flowdir_file) as src:
        extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]

    upstream_runoff_map = funmixer.get_upstream_concentration_map(
        area_dict, runoff_rates
    )

    # --------------------------------
    # Phosphate flux calculation
    # --------------------------------
    conversion_factor = 1e-6 * 4  # convert from ppm kg /m2/ 6 hours to kg/m2/day

    upstream_flux_mean_map = upstream_runoff_map * upstream_mean_map * conversion_factor

    upstream_flux_std_map = upstream_runoff_map * upstream_std_map * conversion_factor

    print("Flux min:", np.nanmin(upstream_flux_mean_map))
    print("Flux max:", np.nanmax(upstream_flux_mean_map))

    # Save upstream_flux_std_map as a GeoTIFF using the same profile as flowdir_file
    with rasterio.open(flowdir_file) as src:
        profile = src.profile.copy()

        # Update profile for the new data
        profile.update(
            dtype=upstream_flux_std_map.dtype,
            count=1,  # Single band
            compress="lzw",  # Optional compression
            # Change from .nc to .tif
            driver="GTiff",
        )
        # Write the GeoTIFF for mean
        output_tif = f"Modelling/DataNEW/phosphate_flux_mean_{season}.tif"
        with rasterio.open(output_tif, "w", **profile) as dst:
            dst.write(upstream_flux_mean_map, 1)

        print(f"Saved GeoTIFF of mean flux → {output_tif}")

        # Write the GeoTIFF for std
        output_tif = f"Modelling/DataNEW/phosphate_flux_std_{season}.tif"
        with rasterio.open(output_tif, "w", **profile) as dst:
            dst.write(upstream_flux_std_map, 1)

        print(f"Saved GeoTIFF of std flux → {output_tif}")

    # -------------------------------
    # Plots
    # -------------------------------

    # -------------------------------
    # River network (once)
    # -------------------------------
    rx, ry, rs = build_river_network(
        flowdir_file,
        drainage_area_threshold,
    )

    plot_upstream_with_rivers(
        upstream_flux_mean_map,
        extent,
        rx,
        ry,
        rs,
        title=title_string,
        cmap="viridis",
        vmin=1e-7,
        vmax=1e-5,
        ax=None,
    )
    # Save this in Figures/ with an appropriate name
    plt.savefig(f"Modelling/Figures/phosphate_flux_mean_{season}.png", dpi=300)
    plt.show()
    plot_upstream_with_rivers(
        upstream_flux_std_map,
        extent,
        rx,
        ry,
        rs,
        title=title_string,
        cmap="magma",
        vmin=1e-8,
        vmax=1e-6,
        ax=None,
    )
    # Save this in Figures/ with an appropriate name
    plt.savefig(f"Modelling/Figures/phosphate_flux_std_{season}.png", dpi=300)
    plt.show()

    return (upstream_flux_mean_map, upstream_flux_std_map, extent, rx, ry, rs)


# ---------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------

flowdir_file = "Data_PP/Original/welsh_d8.nc"

# Check if this directory exist and create if not
output_base = "Modelling/DataNEW"
os.makedirs(output_base, exist_ok=True)

seasons = {
    "Autumn_2023": {
        "unmix": "Modelling/Phosphate_averages/Autumn 2023.csv",
        "runoff": "Runoff/Autumn_2023/runoff_autumn_2023_mean_by_subcatchment.csv",
    },
    "Winter_2023_2024": {
        "unmix": "Modelling/Phosphate_averages/Winter 2023_2024.csv",
        "runoff": "Runoff/Winter_2023_2024/runoff_winter_2023_2024_mean_by_subcatchment.csv",
    },
    "Spring_2024": {
        "unmix": "Modelling/Phosphate_averages/Spring 2024.csv",
        "runoff": "Runoff/Spring_2024/runoff_spring_2024_mean_by_subcatchment.csv",
    },
    "Summer_2024": {
        "unmix": "Modelling/Phosphate_averages/Summer 2024.csv",
        "runoff": "Runoff/Summer_2024/runoff_summer_2024_mean_by_subcatchment.csv",
    },
    "Autumn_2024": {
        "unmix": "Modelling/Phosphate_averages/Autumn 2024.csv",
        "runoff": "Runoff/Autumn_2024/runoff_autumn_2024_mean_by_subcatchment.csv",
    },
    "Winter_2024_2025": {
        "unmix": "Modelling/Phosphate_averages/Winter 2024_2025.csv",
        "runoff": "Runoff/Winter_2024_2025/runoff_winter_2024_2025_mean_by_subcatchment.csv",
    },
    "Spring_2025": {
        "unmix": "Modelling/Phosphate_averages/Spring 2025.csv",
        "runoff": "Runoff/Spring_2025/runoff_spring_2025_mean_by_subcatchment.csv",
    },
    "Summer_2025": {
        "unmix": "Modelling/Phosphate_averages/Summer 2025.csv",
        "runoff": "Runoff/Summer_2025/runoff_summer_2025_mean_by_subcatchment.csv",
    },
}

mean_maps = {}
std_maps = {}

for season, files in seasons.items():
    print("##############################")
    print(f"Processing season: {season}")

    mean_map, std_map, extent, rx, ry, rs = run_funmixer_montecarlo(
        unmix_csv=files["unmix"],
        runoff_csv=files["runoff"],
        flowdir_file=flowdir_file,
        element_col="phosphate_mean",
        season=season,
        output_base=output_base,
    )

    mean_maps[season] = mean_map
    std_maps[season] = std_map

all_mean_values = np.concatenate(
    [m[(~np.isnan(m)) & (m > 0)] for m in mean_maps.values()]
)
all_std_values = np.concatenate(
    [m[(~np.isnan(m)) & (m > 0)] for m in std_maps.values()]
)

mean_vmin = max(np.percentile(all_mean_values, 2), 1e-12)
mean_vmax = np.percentile(all_mean_values, 98)

std_vmin = max(np.percentile(all_std_values, 2), 1e-12)
std_vmax = np.percentile(all_std_values, 98)

fig, axes = plt.subplots(2, 4, figsize=(18, 9), subplot_kw={"projection": ccrs.OSGB()})

axes = axes.flatten()

for ax, (season, mean_map) in zip(axes, mean_maps.items()):
    title = format_season_title(season)

    im = plot_upstream_with_rivers(
        mean_map,
        extent,
        rx,
        ry,
        rs,
        title=title,
        cmap="viridis",
        vmin=mean_vmin,
        vmax=mean_vmax,
        ax=ax,
    )

# ---- Adjust bottom space to fit x-ticks & colorbar ----
fig.subplots_adjust(
    bottom=0.15, top=0.92, left=0.05, right=0.97, hspace=0.3, wspace=0.2
)

cbar = fig.colorbar(im, ax=axes, orientation="horizontal", fraction=0.05, pad=0.1)
cbar.set_label("Phosphate flux (kg/m²/day)")

fig.suptitle("Mean upstream phosphate flux", fontsize=16)
plt.savefig(os.path.join(output_base, "panel_mean_flux.png"), dpi=300)
plt.close()

fig, axes = plt.subplots(2, 4, figsize=(18, 9), subplot_kw={"projection": ccrs.OSGB()})

axes = axes.flatten()

for ax, (season, std_map) in zip(axes, std_maps.items()):
    title = format_season_title(season)

    im = plot_upstream_with_rivers(
        std_map,
        extent,
        rx,
        ry,
        rs,
        title=title,
        cmap="magma",
        vmin=std_vmin,
        vmax=std_vmax,
        ax=ax,
    )

fig.subplots_adjust(
    bottom=0.15, top=0.92, left=0.05, right=0.97, hspace=0.3, wspace=0.2
)
cbar = fig.colorbar(im, ax=axes, orientation="horizontal", fraction=0.05, pad=0.1)
cbar.set_label("Flux uncertainty (kg/m²/day)")

fig.suptitle("Upstream flux uncertainty", fontsize=16)

plt.savefig(os.path.join(output_base, "panel_flux_uncertainty.png"), dpi=300)
plt.close()
