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
from osgeo import gdal
import funmixer.flow_acc_cfuncs as cf
from scipy.stats import gmean as geo_mean

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
        zorder=4,
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

    element_data = funmixer.get_element_obs(element_col, obs)

    #problem = funmixer.SampleNetworkUnmixer(sample_network)
    # -------------------------------
    # Define rate constants
    # -------------------------------
    PHOSPHATE_DECAY_K = 1.45e-6  # m^-1

    rate_constants = {
        site_id: PHOSPHATE_DECAY_K
        for site_id in element_data.keys()
    }

    print("Decay constant range:",
      min(rate_constants.values()),
      max(rate_constants.values()))
      
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
    )

    # Observed vs predicted downstream
    plt.figure(figsize=(5, 5))
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
    plt.savefig(
        "/Users/mahikabhandari/Desktop/Earth Science Year 4/"
        "MSci Project/Analysis/Modelling/Figures/predicted_vs_observed.png",
        dpi=300,
    )
    plt.show()

    def load_runoff_export_rates(
        runoff_csv: str,
        site_col: str = "Site_IDs",
        runoff_col: str = "Mean_Runoff",
        normalise: bool = True,
        min_value: float = 1e-6,
    ) -> dict[str, float]:
        """
        Load runoff CSV (already merged by snapped coordinates)
        and return export rates for funmixer.
        """
        df = pd.read_csv(runoff_csv)

        export_rates = {
            str(row[site_col]): max(float(row[runoff_col]), min_value)
            for _, row in df.iterrows()
        }

        if normalise:
            gm = geo_mean(list(export_rates.values()))
            export_rates = {k: v / gm for k, v in export_rates.items()}

        return export_rates

    runoff_rates = load_runoff_export_rates(
        runoff_csv="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Spring_2024/runoff_spring2024_mean_by_subcatchment.csv"
    )

    missing = set(element_data.keys()) - set(runoff_rates.keys())
    extra = set(runoff_rates.keys()) - set(element_data.keys())

    print("Missing runoff for sites:", missing)
    print("Extra runoff sites:", extra)


    # -------------------------------
    # Monte Carlo solve
    # -------------------------------

    downstream_mc, upstream_mc = problem.solve_montecarlo(
    element_data,
    relative_error=relative_error,
    num_repeats=num_repeats,
    export_rates=runoff_rates,
    #rate_constants=rate_constants,
    regularization_strength=regularization_strength,
    )

    downstream_mean = {k: np.mean(v) for k, v in downstream_mc.items()}
    downstream_std  = {k: np.std(v)  for k, v in downstream_mc.items()}

    upstream_mean = {k: np.mean(v) for k, v in upstream_mc.items()}
    upstream_std  = {k: np.std(v)  for k, v in upstream_mc.items()}

    print("Downstream mean (with decay) min/max:",
      min(downstream_mean.values()),
      max(downstream_mean.values()))

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
        
    upstream_runoff_map = funmixer.get_upstream_concentration_map(
        area_dict, runoff_rates
    ) 

    # --------------------------------
    # Phosphate flux calculation
    # --------------------------------
    CONVERSION_FACTOR = 1e-6 * 4

    upstream_flux_mean_map = (
        upstream_runoff_map * upstream_mean_map * CONVERSION_FACTOR
    )

    upstream_flux_std_map = (
        upstream_runoff_map * upstream_std_map * CONVERSION_FACTOR
    )

    print("Flux min:", np.nanmin(upstream_flux_mean_map))
    print("Flux max:", np.nanmax(upstream_flux_mean_map))


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
    upstream_flux_mean_map,
    extent,
    rx, ry, rs,
    title="Upstream Mean Phosphate Flux",
    cbar_label="Phosphate flux (kg/m²/day)",
    cmap="viridis",
    output_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/"
    "MSci Project/Analysis/Modelling/Figures/upstream_flux_mean_montecarlo.png",
    )

    plot_upstream_with_rivers(
        upstream_flux_std_map,
        extent,
        rx, ry, rs,
        title="Upstream Phosphate Flux Uncertainty",
        cbar_label="Flux uncertainty (kg/m²/day)",
        cmap="magma",
        output_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/"
        "MSci Project/Analysis/Modelling/Figures/upstream_flux_uncertainty_montecarlo.png",
    )

# ---------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------
run_funmixer_montecarlo(
    unmix_csv="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Modelling/phosphate_averages.csv",
    flowdir_file="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Original/welsh_d8.nc",
    element_col="Phosphate_mean",
)

