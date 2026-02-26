"""Reprojects EFAS runoff grids onto same coordinate system as the flow-direction grid, and calculates average seasonal runoff."""

import rioxarray as rxr
import xarray as xr
from rasterio.enums import Resampling


def runoff_seasonal_mean_to_d8(
    runoff_nc_path,
    d8_nc_path,
    output_nc_path,
    season_months,
    runoff_var="rowe",
    runoff_crs="EPSG:4326",
    d8_crs="EPSG:27700",
):
    """
    Compute seasonal mean runoff and reproject to match Welsh D8 grid.

    Parameters
    ----------
    runoff_nc_path : str
        Path to runoff NetCDF file
    d8_nc_path : str
        Path to welsh_d8.nc
    output_nc_path : str
        Output NetCDF path
    season_months : list
        List of month numbers for the season (e.g., [3, 4, 5] for MAM)
    runoff_var : str
        Runoff variable name
    runoff_crs : str
        CRS of runoff dataset (default: EPSG:4326)
    d8_crs : str
        CRS of D8 dataset (default: EPSG:27700)

    """
    # Open runoff dataset
    ds_runoff = xr.open_dataset(runoff_nc_path)
    runoff = ds_runoff[runoff_var]

    # Rename time coordinate if needed (GloFAS quirk)
    if "valid_time" in runoff.coords:
        runoff = runoff.rename({"valid_time": "time"})

    # Set spatial dimensions and CRS
    runoff = runoff.rio.set_spatial_dims(
        x_dim="longitude", y_dim="latitude", inplace=False
    )

    if runoff.rio.crs is None:
        runoff = runoff.rio.write_crs(runoff_crs, inplace=False)

    print(f"Runoff CRS: {runoff.rio.crs}")

    # Select seasonal months and compute mean
    runoff_seasonal = runoff.sel(time=runoff["time"].dt.month.isin(season_months))
    runoff_seasonal_mean = runoff_seasonal.mean(dim="time", skipna=True)

    # Open D8 dataset (target grid)
    d8_ds = xr.open_dataset(d8_nc_path)
    d8_da = next(iter(d8_ds.data_vars.values()))

    d8_da = d8_da.rio.set_spatial_dims(
        x_dim=d8_da.rio.x_dim, y_dim=d8_da.rio.y_dim, inplace=False
    )

    if d8_da.rio.crs is None:
        d8_da = d8_da.rio.write_crs(d8_crs, inplace=False)

    print(f"D8 CRS: {d8_da.rio.crs}")

    # Reproject to match D8 grid - Could use bilinear resampling as it is continuous
    runoff_on_d8 = runoff_seasonal_mean.rio.reproject_match(
        d8_da, resampling=Resampling.nearest
    )

    # Save output as both NetCDF and GeoTIFF
    runoff_on_d8.to_netcdf(output_nc_path)
    output_tif_path = output_nc_path.replace(".nc", ".tif")
    runoff_on_d8.rio.to_raster(output_tif_path)
    print(f"✅ Saved NetCDF to: {output_nc_path}")
    print(f"✅ Saved GeoTIFF to: {output_tif_path}")

    return runoff_on_d8


def validate_output(output_path):
    """Validate the output file CRS using both xarray and rioxarray."""
    out_1 = xr.open_dataset(output_path)
    out_2 = rxr.open_rasterio(output_path)

    print(f"Dataset: {out_1}")
    print(f"XArray CRS: {out_1.rio.crs}")
    print(f"RioXArray CRS: {out_2.rio.crs}")


# Configuration for all seasons and years
d8_path = "Data_PP/Original/welsh_d8.nc"

seasons_config = [
    # Season, months, years, season_name
    ("Spring", [3, 4, 5], [2024, 2025], "spring"),
    ("Summer", [6, 7, 8], [2024, 2025], "summer"),
    ("Autumn", [9, 10, 11], [2023, 2024], "autumn"),
    ("Winter", [12, 1, 2], ["2023_2024", "2024_2025"], "winter"),
]

# Process all seasons and years
for season, months, years, season_name in seasons_config:
    print(f"\n{'='*50}")
    print(f"Processing {season} datasets")
    print(f"{'='*50}")

    for year in years:
        if season == "Winter":
            input_path = f"Runoff/{season}_{year}/runoff_{season_name}_{year}.nc"
        else:
            input_path = f"Runoff/{season}_{year}/runoff_{season_name}{year}.nc"
        output_path = (
            f"Runoff/{season}_{year}/runoff_{season_name}_{year}_mean_wales.nc"
        )
        print(f"\nProcessing {season} {year}...")

        # Process the data
        runoff_seasonal_mean_to_d8(
            runoff_nc_path=input_path,
            d8_nc_path=d8_path,
            output_nc_path=output_path,
            season_months=months,
        )

        # # Validate the output
        # validate_output(output_path) # Commenting this out as it causes a segmentation fault for me!
