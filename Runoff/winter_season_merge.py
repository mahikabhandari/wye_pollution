import xarray as xr

def merge_winter_nc(dec_nc, jf_nc, output_nc):
    """
    Merge December (year N) and January–February (year N+1) NetCDF files
    into a single winter (DJF) NetCDF file.

    Parameters
    ----------
    dec_nc : str
        Path to December NetCDF file (e.g. Dec 2023)
    jf_nc : str
        Path to January–February NetCDF file (e.g. Jan–Feb 2024)
    output_nc : str
        Path to save merged winter NetCDF file
    """
    ds = xr.open_mfdataset(
        [dec_nc, jf_nc],
        combine="by_coords"
    )

    # Optional: sort by time just to be safe
    ds = ds.sortby("valid_time")

    # Save merged dataset
    ds.to_netcdf(output_nc)

    print(f"Winter dataset saved to:\n{output_nc}")


dec_nc = "/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Winter_2024_2025/runoff_dec_2024.nc"
jf_nc = "/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Winter_2024_2025/runoff_jan_feb_2025.nc"
output_nc = "/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Winter_2024_2025/runoff_winter_2024_2025.nc"

merge_winter_nc(dec_nc, jf_nc, output_nc)

