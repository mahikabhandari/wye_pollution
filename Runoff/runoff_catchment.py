"""Calculate average runoff per subcatchment of River Wye."""
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import numpy.typing as npt
import pandas as pd
from funmixer import get_sample_graph, get_upstream_concentration_map
from funmixer.d8processing import read_geo_file, write_geotiff

# ---------------------- Build sample network and labels ----------------------------

flowdir_file = "Data_PP/Original/welsh_d8.nc"


def get_unique_upstream_areas(
    sample_network: nx.DiGraph,
    labels: npt.NDArray[np.int_],
) -> dict[str, npt.NDArray[np.bool_]]:
    """
    Get unique upstream areas for each node in the sample network.

    Parameters
    ----------
    sample_network : nx.DiGraph
        The sample network graph.
    labels : npt.NDArray[np.int_]
        The array of subcatchment labels.

    Returns
    -------
    dict[str, npt.NDArray[np.bool_]]
        A dictionary mapping node IDs to boolean masks of their upstream areas.

    """
    return {
        node: labels == node_data["data"].label
        for node, node_data in sample_network.nodes(data=True)
    }


def average_runoff_per_subcatchment(
    sample_network: nx.DiGraph,
    labels: npt.NDArray[np.int_],
    runoff_raster_path: str,
    nodata: float | None = None,
) -> dict[str, float]:
    """
    Calculate the average runoff for each subcatchment in the sample network.

    Parameters
    ----------
    sample_network : nx.DiGraph
        The sample network graph.
    labels : npt.NDArray[np.int_]
        The array of subcatchment labels.
    runoff_raster_path : str
        Path to the raster file containing runoff data.
    nodata : float | None
        No-data value to ignore when calculating averages.

    Returns
    -------
    dict[str, float]
        A dictionary mapping subcatchment IDs to their average runoff values.

    """
    runoff, ds = read_geo_file(runoff_raster_path)
    # Flip runoff upside down to match the labels tif
    runoff = np.flipud(runoff)

    # plt.imshow(runoff)
    # plt.colorbar(label="Runoff")
    # plt.show()

    # plt.imshow(labels)
    # plt.colorbar(label="Subcatchment Labels")
    # plt.show()

    print("Runoff shape:", runoff.shape)
    print("Labels shape:", labels.shape)

    upstream_masks = get_unique_upstream_areas(sample_network, labels)

    mean_runoff: dict[str, float] = {}

    for sample_id, mask in upstream_masks.items():
        values = np.array(runoff[mask], copy=True)

        if nodata is not None:
            values = values[values != nodata]

        mean_runoff[sample_id] = (
            float(np.nanmean(values)) if values.size > 0 else np.nan
        )

    mean_runoff_map = get_upstream_concentration_map(upstream_masks, mean_runoff)

    # Create a 3 part multipanel plot that shows the runoff, the labels and the average so we can sense check everything!
    fig, axs = plt.subplots(1, 3, figsize=(12, 3))
    runoff_min, runoff_max = np.nanmin(runoff), np.nanmax(runoff)
    axs[0].imshow(runoff, cmap="Blues", vmin=runoff_min, vmax=runoff_max)
    # Add colorbar

    cb = plt.colorbar(
        axs[0].imshow(runoff, vmin=runoff_min, vmax=runoff_max, cmap="Blues"), ax=axs[0]
    )
    axs[0].set_xlim(2000, 4500)
    axs[0].set_ylim(4800, 2000)
    cb.set_label("Runoff")
    axs[0].set_title("Runoff")
    # Add a grid
    axs[0].grid(True)

    axs[1].imshow(labels)
    axs[1].set_xlim(2000, 4500)
    axs[1].set_ylim(4800, 2000)
    axs[1].set_title("Subcatchment Labels")
    axs[1].grid(True)

    axs[2].imshow(mean_runoff_map, cmap="Blues", vmin=runoff_min, vmax=runoff_max)
    axs[2].set_xlim(2000, 4500)
    axs[2].set_ylim(4800, 2000)
    cb = plt.colorbar(
        axs[2].imshow(mean_runoff_map, cmap="Blues", vmin=runoff_min, vmax=runoff_max),
        ax=axs[2],
    )
    cb.set_label("Average Runoff")
    axs[2].set_title("Average Runoff")
    axs[2].grid(True)

    plt.tight_layout()
    plt.show()

    return mean_runoff, (mean_runoff_map, ds)


def process_seasonal_runoff(season_name: str, year: str):
    """
    Process runoff data for a given season and year.

    Parameters
    ----------
    season_name : str
        Name of the season (e.g., 'spring', 'summer')
    year : str
        Year identifier (e.g., '2024', '2023_2024')
    base_path : str
        Base path for data files

    """
    print(f"\n{'='*60}")
    print(f"Processing {season_name.title()} {year} runoff data")
    print(f"{'='*60}")

    # Construct file paths
    runoff_raster_path = f"Runoff/{season_name.title()}_{year}/runoff_{season_name}_{year}_mean_wales.tif"
    out_csv = f"Runoff/{season_name.title()}_{year}/runoff_{season_name}_{year}_mean_by_subcatchment.csv"
    input_phosphate_data = (
        f"Modelling/Phosphate_averages/{season_name.title()} {year}.csv"
    )

    sample_network, labels = get_sample_graph(
        flowdir_file,
        input_phosphate_data,
    )

    n_nodes = sample_network.number_of_nodes()
    print(f"Sample network loaded with {n_nodes} nodes.")

    # Calculate mean runoff per subcatchment
    mean_runoff, mean_runoff_map_ds = average_runoff_per_subcatchment(
        sample_network=sample_network,
        labels=labels,
        runoff_raster_path=runoff_raster_path,
        nodata=-9999,
    )

    # Save the mean_runoff_map as a tif in the output directory
    output_tif_path = f"Runoff/{season_name.title()}_{year}/mean_runoff_by_subcatchment_{season_name}_{year}.tif"
    write_geotiff(output_tif_path, mean_runoff_map_ds[0], mean_runoff_map_ds[1])
    print(f"✅ Saved mean runoff map to:\n{output_tif_path}")

    # Create DataFrame and save to CSV
    runoff_df = pd.DataFrame(
        {
            "Site_IDs": list(mean_runoff.keys()),
            "Mean_Runoff": list(mean_runoff.values()),
        }
    )

    runoff_df.to_csv(out_csv, index=False)
    print(
        f"✅ Saved {season_name.title()} {year} mean runoff per subcatchment to:\n{out_csv}"
    )


seasons_config = [
    ("spring", "2024"),
    ("spring", "2025"),
    ("summer", "2024"),
    ("summer", "2025"),
    ("autumn", "2023"),
    ("autumn", "2024"),
    ("winter", "2023_2024"),
    ("winter", "2024_2025"),
]

# Process all seasons and years
for season, year in seasons_config:
    process_seasonal_runoff(season, year)

print(f"\n{'='*60}")
print("All seasonal runoff processing complete!")
print(f"{'='*60}")
