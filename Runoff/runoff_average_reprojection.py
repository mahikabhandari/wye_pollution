import xarray as xr
import rioxarray as rxr
from rasterio.enums import Resampling
import rioxarray as rxr

"""
Spring 2024 and 2024 runoff dataset
"""
ds = xr.open_dataset("/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Spring_2024/runoff_spring2024.nc")
print(ds)

def runoff_MAM_mean_to_d8(
    runoff_nc_path,
    d8_nc_path,
    output_nc_path,
    runoff_var="runoff_water_equivalent",
    runoff_crs="EPSG:4326",
    d8_crs="EPSG:27700",
):
    """
    Compute March–May mean runoff 2024 and reproject to match Welsh D8 grid.

    Parameters
    ----------
    runoff_nc_path : str
        Path to GloFAS NetCDF (or wildcard for multiple files)
    d8_nc_path : str
        Path to welsh_d8.nc
    output_nc_path : str
        Output NetCDF path
    runoff_var : str
        Runoff variable name
    runoff_crs : str
        CRS of runoff dataset (default: EPSG:4326)
    d8_crs : str
        CRS of D8 dataset (default: EPSG:27700)
    """

    # ------------------------------------------------------------------
    # Open runoff dataset
    # ------------------------------------------------------------------
    ds_runoff = xr.open_dataset(runoff_nc_path)
    runoff = ds_runoff[runoff_var]

    # Rename time coordinate if needed (GloFAS quirk)
    if "valid_time" in runoff.coords:
        runoff = runoff.rename({"valid_time": "time"})

    # Explicitly define spatial dimensions (VERY IMPORTANT)
    runoff = runoff.rio.set_spatial_dims(
        x_dim="longitude",
        y_dim="latitude",
        inplace=False,
    )

    # Attach CRS ONLY if missing
    if runoff.rio.crs is None:
        runoff = runoff.rio.write_crs(runoff_crs, inplace=False)

    print(f"Runoff CRS: {runoff.rio.crs}")

    # ------------------------------------------------------------------
    # Select March–May and compute seasonal mean
    # ------------------------------------------------------------------
    runoff_mam = runoff.sel(
        time=runoff["time"].dt.month.isin([3, 4, 5])
    )

    runoff_mam_mean = runoff_mam.mean(dim="time", skipna=True)

    # ------------------------------------------------------------------
    # Open D8 dataset (this defines the target grid)
    # ------------------------------------------------------------------
    d8_ds = xr.open_dataset(d8_nc_path)

    # Assume only one data variable in D8 file
    d8_da = next(iter(d8_ds.data_vars.values()))

    # Explicit spatial dims for D8
    print(f"D8 original dims: {d8_da.rio.x_dim}, {d8_da.rio.y_dim}")
    d8_da = d8_da.rio.set_spatial_dims(
        x_dim=d8_da.rio.x_dim,
        y_dim=d8_da.rio.y_dim,
        inplace=False,
    )

    # Attach CRS ONLY if missing
    if d8_da.rio.crs is None:
        d8_da = d8_da.rio.write_crs(d8_crs, inplace=False)

    print(f"D8 CRS: {d8_da.rio.crs}")

    # ------------------------------------------------------------------
    # Reproject and resample to EXACT D8 grid
    # ------------------------------------------------------------------
    runoff_on_d8 = runoff_mam_mean.rio.reproject_match(
        d8_da,
        #resampling=Resampling.bilinear  # appropriate for continuous data
        resampling=Resampling.nearest  # appropriate for categorical data
    )

    # ------------------------------------------------------------------
    # Save output
    # ------------------------------------------------------------------
    runoff_on_d8.to_netcdf(output_nc_path)

    print("✅ Runoff successfully reprojected to Welsh D8 grid")
    print(f"📁 Saved to: {output_nc_path}")

    return runoff_on_d8

runoff_MAM_mean_to_d8(
    runoff_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Spring_2024/runoff_spring2024.nc",
    d8_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Original/welsh_d8.nc",
    output_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Spring_2024/runoff_spring2024_mean_wales.nc",
    runoff_var="rowe"
)

out_1 = xr.open_dataset("/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Spring_2024/runoff_spring2024_mean_wales.nc")

print(out_1)
print(out_1.rio.crs)

out_2 = rxr.open_rasterio(
    "/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Spring_2024/runoff_spring2024_mean_wales.nc"
)

print(out_2.rio.crs)


runoff_MAM_mean_to_d8(
    runoff_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Spring_2025/runoff_spring2025.nc",
    d8_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Original/welsh_d8.nc",
    output_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Spring_2025/runoff_spring2025_mean_wales.nc",
    runoff_var="rowe"
)

out_1 = xr.open_dataset("/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Spring_2025/runoff_spring2025_mean_wales.nc")

print(out_1)
print(out_1.rio.crs)

out_2 = rxr.open_rasterio(
    "/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Spring_2025/runoff_spring2025_mean_wales.nc"
)

print(out_2.rio.crs)


"""
Summer 2024 and 2025 runoff dataset
"""

def runoff_JJA_mean_to_d8(
    runoff_nc_path,
    d8_nc_path,
    output_nc_path,
    runoff_var="runoff_water_equivalent",
    runoff_crs="EPSG:4326",
    d8_crs="EPSG:27700",
):
    """
    Compute June–August mean runoff and reproject to match Welsh D8 grid.

    Parameters
    ----------
    runoff_nc_path : str
        Path to GloFAS NetCDF (or wildcard for multiple files)
    d8_nc_path : str
        Path to welsh_d8.nc
    output_nc_path : str
        Output NetCDF path
    runoff_var : str
        Runoff variable name
    runoff_crs : str
        CRS of runoff dataset (default: EPSG:4326)
    d8_crs : str
        CRS of D8 dataset (default: EPSG:27700)
    """

    # ------------------------------------------------------------------
    # Open runoff dataset
    # ------------------------------------------------------------------
    ds_runoff = xr.open_dataset(runoff_nc_path)
    runoff = ds_runoff[runoff_var]

    # Rename time coordinate if needed (GloFAS quirk)
    if "valid_time" in runoff.coords:
        runoff = runoff.rename({"valid_time": "time"})

    # Explicitly define spatial dimensions (VERY IMPORTANT)
    runoff = runoff.rio.set_spatial_dims(
        x_dim="longitude",
        y_dim="latitude",
        inplace=False,
    )

    # Attach CRS ONLY if missing
    if runoff.rio.crs is None:
        runoff = runoff.rio.write_crs(runoff_crs, inplace=False)

    print(f"Runoff CRS: {runoff.rio.crs}")

    # ------------------------------------------------------------------
    # Select March–May and compute seasonal mean
    # ------------------------------------------------------------------
    runoff_jja = runoff.sel(
        time=runoff["time"].dt.month.isin([6, 7, 8])
    )

    runoff_jja_mean = runoff_jja.mean(dim="time", skipna=True)

    # ------------------------------------------------------------------
    # Open D8 dataset (this defines the target grid)
    # ------------------------------------------------------------------
    d8_ds = xr.open_dataset(d8_nc_path)

    # Assume only one data variable in D8 file
    d8_da = next(iter(d8_ds.data_vars.values()))

    # Explicit spatial dims for D8
    print(f"D8 original dims: {d8_da.rio.x_dim}, {d8_da.rio.y_dim}")
    d8_da = d8_da.rio.set_spatial_dims(
        x_dim=d8_da.rio.x_dim,
        y_dim=d8_da.rio.y_dim,
        inplace=False,
    )

    # Attach CRS ONLY if missing
    if d8_da.rio.crs is None:
        d8_da = d8_da.rio.write_crs(d8_crs, inplace=False)

    print(f"D8 CRS: {d8_da.rio.crs}")

    # ------------------------------------------------------------------
    # Reproject and resample to EXACT D8 grid
    # ------------------------------------------------------------------
    runoff_on_d8 = runoff_jja_mean.rio.reproject_match(
        d8_da,
        #resampling=Resampling.bilinear  # appropriate for continuous data
        resampling=Resampling.nearest  # appropriate for categorical data
    )

    # ------------------------------------------------------------------
    # Save output
    # ------------------------------------------------------------------
    runoff_on_d8.to_netcdf(output_nc_path)

    print("✅ Runoff successfully reprojected to Welsh D8 grid")
    print(f"📁 Saved to: {output_nc_path}")

    return runoff_on_d8

runoff_JJA_mean_to_d8(
    runoff_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Summer_2024/runoff_summer2024.nc",
    d8_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Original/welsh_d8.nc",
    output_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Summer_2024/runoff_summer2024_mean_wales.nc",
    runoff_var="rowe"
)

out_1 = xr.open_dataset("/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Summer_2024/runoff_summer2024_mean_wales.nc")

print(out_1)
print(out_1.rio.crs)

out_2 = rxr.open_rasterio(
    "/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Summer_2024/runoff_summer2024_mean_wales.nc"
)

print(out_2.rio.crs)

runoff_JJA_mean_to_d8(
    runoff_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Summer_2025/runoff_summer2025.nc",
    d8_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Original/welsh_d8.nc",
    output_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Summer_2025/runoff_summer2025_mean_wales.nc",
    runoff_var="rowe"
)

out_1 = xr.open_dataset("/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Summer_2025/runoff_summer2025_mean_wales.nc")

print(out_1)
print(out_1.rio.crs)

out_2 = rxr.open_rasterio(
    "/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Summer_2025/runoff_summer2025_mean_wales.nc"
)

print(out_2.rio.crs)

"""
Autumn 2023 and 2024 runoff dataset
"""

def runoff_SON_mean_to_d8(
    runoff_nc_path,
    d8_nc_path,
    output_nc_path,
    runoff_var="runoff_water_equivalent",
    runoff_crs="EPSG:4326",
    d8_crs="EPSG:27700",
):
    """
    Compute June–August mean runoff and reproject to match Welsh D8 grid.

    Parameters
    ----------
    runoff_nc_path : str
        Path to GloFAS NetCDF (or wildcard for multiple files)
    d8_nc_path : str
        Path to welsh_d8.nc
    output_nc_path : str
        Output NetCDF path
    runoff_var : str
        Runoff variable name
    runoff_crs : str
        CRS of runoff dataset (default: EPSG:4326)
    d8_crs : str
        CRS of D8 dataset (default: EPSG:27700)
    """

    # ------------------------------------------------------------------
    # Open runoff dataset
    # ------------------------------------------------------------------
    ds_runoff = xr.open_dataset(runoff_nc_path)
    runoff = ds_runoff[runoff_var]

    # Rename time coordinate if needed (GloFAS quirk)
    if "valid_time" in runoff.coords:
        runoff = runoff.rename({"valid_time": "time"})

    # Explicitly define spatial dimensions (VERY IMPORTANT)
    runoff = runoff.rio.set_spatial_dims(
        x_dim="longitude",
        y_dim="latitude",
        inplace=False,
    )

    # Attach CRS ONLY if missing
    if runoff.rio.crs is None:
        runoff = runoff.rio.write_crs(runoff_crs, inplace=False)

    print(f"Runoff CRS: {runoff.rio.crs}")

    # ------------------------------------------------------------------
    # Select September–November and compute seasonal mean
    # ------------------------------------------------------------------
    runoff_son_mean = runoff.sel(
        time=runoff["time"].dt.month.isin([9, 10, 11])
    )

    runoff_son_mean = runoff_son_mean.mean(dim="time", skipna=True)

    # ------------------------------------------------------------------
    # Open D8 dataset (this defines the target grid)
    # ------------------------------------------------------------------
    d8_ds = xr.open_dataset(d8_nc_path)

    # Assume only one data variable in D8 file
    d8_da = next(iter(d8_ds.data_vars.values()))

    # Explicit spatial dims for D8
    print(f"D8 original dims: {d8_da.rio.x_dim}, {d8_da.rio.y_dim}")
    d8_da = d8_da.rio.set_spatial_dims(
        x_dim=d8_da.rio.x_dim,
        y_dim=d8_da.rio.y_dim,
        inplace=False,
    )

    # Attach CRS ONLY if missing
    if d8_da.rio.crs is None:
        d8_da = d8_da.rio.write_crs(d8_crs, inplace=False)

    print(f"D8 CRS: {d8_da.rio.crs}")

    # ------------------------------------------------------------------
    # Reproject and resample to EXACT D8 grid
    # ------------------------------------------------------------------
    runoff_on_d8 = runoff_son_mean.rio.reproject_match(
        d8_da,
        #resampling=Resampling.bilinear  # appropriate for continuous data
        resampling=Resampling.nearest  # appropriate for categorical data
    )

    # ------------------------------------------------------------------
    # Save output
    # ------------------------------------------------------------------
    runoff_on_d8.to_netcdf(output_nc_path)

    print("✅ Runoff successfully reprojected to Welsh D8 grid")
    print(f"📁 Saved to: {output_nc_path}")

    return runoff_on_d8


runoff_SON_mean_to_d8(
    runoff_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Autumn_2023/runoff_autumn2023.nc",
    d8_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Original/welsh_d8.nc",
    output_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Autumn_2023/runoff_autumn2023_mean_wales.nc",
    runoff_var="rowe"
)

out_1 = xr.open_dataset("/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Autumn_2023/runoff_autumn2023_mean_wales.nc")

print(out_1)
print(out_1.rio.crs)

out_2 = rxr.open_rasterio(
    "/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Autumn_2023/runoff_autumn2023_mean_wales.nc"
)

print(out_2.rio.crs)

runoff_SON_mean_to_d8(
    runoff_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Autumn_2024/runoff_autumn2024.nc",
    d8_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Original/welsh_d8.nc",
    output_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Autumn_2024/runoff_autumn2024_mean_wales.nc",
    runoff_var="rowe"
)

out_1 = xr.open_dataset("/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Autumn_2024/runoff_autumn2024_mean_wales.nc")

print(out_1)
print(out_1.rio.crs)

out_2 = rxr.open_rasterio(
    "/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Autumn_2024/runoff_autumn2024_mean_wales.nc"
)

print(out_2.rio.crs)


"""
Winter 2023/2024 and 2024/2025 runoff dataset
"""

def runoff_DJF_mean_to_d8(
    runoff_nc_path,
    d8_nc_path,
    output_nc_path,
    runoff_var="runoff_water_equivalent",
    runoff_crs="EPSG:4326",
    d8_crs="EPSG:27700",
):
    """
    Compute December-February mean runoff and reproject to match Welsh D8 grid.

    Parameters
    ----------
    runoff_nc_path : str
        Path to GloFAS NetCDF (or wildcard for multiple files)
    d8_nc_path : str
        Path to welsh_d8.nc
    output_nc_path : str
        Output NetCDF path
    runoff_var : str
        Runoff variable name
    runoff_crs : str
        CRS of runoff dataset (default: EPSG:4326)
    d8_crs : str
        CRS of D8 dataset (default: EPSG:27700)
    """

    # ------------------------------------------------------------------
    # Open runoff dataset
    # ------------------------------------------------------------------
    ds_runoff = xr.open_dataset(runoff_nc_path)
    runoff = ds_runoff[runoff_var]

    # Rename time coordinate if needed (GloFAS quirk)
    if "valid_time" in runoff.coords:
        runoff = runoff.rename({"valid_time": "time"})

    # Explicitly define spatial dimensions (VERY IMPORTANT)
    runoff = runoff.rio.set_spatial_dims(
        x_dim="longitude",
        y_dim="latitude",
        inplace=False,
    )

    # Attach CRS ONLY if missing
    if runoff.rio.crs is None:
        runoff = runoff.rio.write_crs(runoff_crs, inplace=False)

    print(f"Runoff CRS: {runoff.rio.crs}")

    # ------------------------------------------------------------------
    # Select September–November and compute seasonal mean
    # ------------------------------------------------------------------
    runoff_son_mean = runoff.sel(
        time=runoff["time"].dt.month.isin([12, 1, 2])
    )

    runoff_son_mean = runoff_son_mean.mean(dim="time", skipna=True)

    # ------------------------------------------------------------------
    # Open D8 dataset (this defines the target grid)
    # ------------------------------------------------------------------
    d8_ds = xr.open_dataset(d8_nc_path)

    # Assume only one data variable in D8 file
    d8_da = next(iter(d8_ds.data_vars.values()))

    # Explicit spatial dims for D8
    print(f"D8 original dims: {d8_da.rio.x_dim}, {d8_da.rio.y_dim}")
    d8_da = d8_da.rio.set_spatial_dims(
        x_dim=d8_da.rio.x_dim,
        y_dim=d8_da.rio.y_dim,
        inplace=False,
    )

    # Attach CRS ONLY if missing
    if d8_da.rio.crs is None:
        d8_da = d8_da.rio.write_crs(d8_crs, inplace=False)

    print(f"D8 CRS: {d8_da.rio.crs}")

    # ------------------------------------------------------------------
    # Reproject and resample to EXACT D8 grid
    # ------------------------------------------------------------------
    runoff_on_d8 = runoff_son_mean.rio.reproject_match(
        d8_da,
        #resampling=Resampling.bilinear  # appropriate for continuous data
        resampling=Resampling.nearest  # appropriate for categorical data
    )

    # ------------------------------------------------------------------
    # Save output
    # ------------------------------------------------------------------
    runoff_on_d8.to_netcdf(output_nc_path)

    print("✅ Runoff successfully reprojected to Welsh D8 grid")
    print(f"📁 Saved to: {output_nc_path}")

    return runoff_on_d8


runoff_DJF_mean_to_d8(
    runoff_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Winter_2023_2024/runoff_winter_2023_2024.nc",
    d8_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Original/welsh_d8.nc",
    output_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Winter_2023_2024/runoff_winter_2023_2024_mean_wales.nc",
    runoff_var="rowe"
)

out_1 = xr.open_dataset("/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Winter_2023_2024/runoff_winter_2023_2024_mean_wales.nc")

print(out_1)
print(out_1.rio.crs)

out_2 = rxr.open_rasterio(
    "/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Winter_2023_2024/runoff_winter_2023_2024_mean_wales.nc"
)

print(out_2.rio.crs)

runoff_DJF_mean_to_d8(
    runoff_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Winter_2024_2025/runoff_winter_2024_2025.nc",
    d8_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Original/welsh_d8.nc",
    output_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Winter_2024_2025/runoff_winter_2024_2025_mean_wales.nc",
    runoff_var="rowe"
)

out_1 = xr.open_dataset("/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Winter_2024_2025/runoff_winter_2024_2025_mean_wales.nc")

print(out_1)
print(out_1.rio.crs)

out_2 = rxr.open_rasterio(
    "/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Winter_2024_2025/runoff_winter_2024_2025_mean_wales.nc"
)

print(out_2.rio.crs)