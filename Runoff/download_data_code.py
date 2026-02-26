"""Downloads EFAS runoff data and over the Wye region and saves it into specified directories for each 3 month season."""

import os
import zipfile

import cdsapi
import xarray as xr


def extract_nc_from_zip(zip_path, output_nc_path):
    """
    Extract the first .nc file from a zip archive and saves it to a specified path.

    Parameters
    ----------
    zip_path : str
        Path to the .zip file
    output_nc_path : str
        Desired path for the extracted .nc file

    """
    with zipfile.ZipFile(zip_path, "r") as zip_ref:
        # Find .nc files in the zip
        nc_files = [f for f in zip_ref.namelist() if f.endswith(".nc")]

        if not nc_files:
            raise FileNotFoundError("No .nc file found in the zip archive.")

        # Extract the first .nc file
        nc_file = nc_files[0]
        zip_ref.extract(nc_file, os.path.dirname(output_nc_path))

        # Rename/move to desired output path
        extracted_path = os.path.join(os.path.dirname(output_nc_path), nc_file)
        os.rename(extracted_path, output_nc_path)

    print(f"Saved NetCDF file to: {output_nc_path}")


def merge_winter_nc(dec_nc: str, jf_nc: str, output_nc: str):
    """
    Merge December (year N) and January-February (year N+1) NetCDF files into a single winter (DJF) NetCDF file.

    Parameters
    ----------
    dec_nc : str
        Path to December NetCDF file (e.g. Dec 2023)
    jf_nc : str
        Path to January-February NetCDF file (e.g. Jan-Feb 2024)
    output_nc : str
        Path to save merged winter NetCDF file

    """
    ds = xr.open_mfdataset([dec_nc, jf_nc], combine="by_coords")

    # Optional: sort by time just to be safe
    ds = ds.sortby("valid_time")

    # Save merged dataset
    ds.to_netcdf(output_nc)

    print(f"Winter dataset saved to:\n{output_nc}")


dataset = "efas-historical"
base_request = {
    "system_version": ["version_5_0"],
    "variable": ["runoff_water_equivalent"],
    "model_levels": "surface_level",
    "hyear": None,  # These will be populated later!
    "hmonth": None,  # These will be populated later!
    "hday": [
        "01",
        "02",
        "03",
        "04",
        "05",
        "06",
        "07",
        "08",
        "09",
        "10",
        "11",
        "12",
        "13",
        "14",
        "15",
        "16",
        "17",
        "18",
        "19",
        "20",
        "21",
        "22",
        "23",
        "24",
        "25",
        "26",
        "27",
        "28",
        "29",
        "30",
        "31",
    ],
    "time": ["00:00", "06:00", "12:00", "18:00"],
    "data_format": "netcdf",
    "download_format": "zip",
    "area": [52.6, -3.9, 51.3, -2.3],
}

# Create a dictionary mapping time period onto output dir
output_dirs = {
    "Autumn_2023": "Runoff/Autumn_2023/",
    "Dec_2023": "Runoff/Winter_2023_2024/",
    "JanFeb_2024": "Runoff/Winter_2023_2024/",
    "Spring_2024": "Runoff/Spring_2024/",
    "Summer_2024": "Runoff/Summer_2024/",
    "Autumn_2024": "Runoff/Autumn_2024/",
    "Dec_2024": "Runoff/Winter_2024_2025/",
    "JanFeb_2025": "Runoff/Winter_2024_2025/",
    "Spring_2025": "Runoff/Spring_2025/",
    "Summer_2025": "Runoff/Summer_2025/",
}

# For each output directory check if it exists and create it if not

for outputdir in output_dirs.values():
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

print("#" * 50)
print("Downloading runoff data, this can take sometime!")

for period, outputdir in output_dirs.items():
    months, year = period.split("_")
    request = base_request.copy()
    if months == "Dec":
        request["hmonth"] = ["12"]
    if months == "JanFeb":
        request["hmonth"] = ["01", "02"]
    if months == "Spring":
        request["hmonth"] = ["03", "04", "05"]
    if months == "Summer":
        request["hmonth"] = ["06", "07", "08"]
    if months == "Autumn":
        request["hmonth"] = ["09", "10", "11"]
    request["hyear"] = [year]

    months_lowercase = months.lower()
    client = cdsapi.Client()
    if months != "Dec" and months != "JanFeb":
        basefilename = f"Runoff/{period}/runoff_{months_lowercase}{year}"
        filename_zip = f"{basefilename}.zip"
        # nc filename is same but changing .zip to .nc
        filename_nc = f"{basefilename}.nc"
    elif months == "Dec":
        basefilename = f"Runoff/Winter_{year}_{int(year) + 1}/runoff_dec_{year}"
        filename_zip = f"{basefilename}.zip"
        filename_nc = f"{basefilename}.nc"
    elif months == "JanFeb":
        basefilename = f"Runoff/Winter_{int(year)-1}_{int(year)}/runoff_jan_feb_{year}"
        filename_zip = f"{basefilename}.zip"
        filename_nc = f"{basefilename}.nc"

    # To test the logic is correct, print which period you are querying, and where it is saving everything to:
    print(
        f"Requesting data for {period} and saving to {filename_zip} and {filename_nc}"
    )
    client.retrieve(dataset, request).download(filename_zip)
    print(f"Downloaded zip file for {period} to {filename_zip}")
    print(f"Extracting NetCDF file for {period} to {filename_nc}")
    extract_nc_from_zip(filename_zip, filename_nc)

print("Finished extracting NetCDF files!")
print("#" * 50)

# Merge the winter seasons which are in different years into one seasonal file
winter_seasons = [("2023", "2024"), ("2024", "2025")]

print("Merging the winter seasons...")

for dec_year, jf_year in winter_seasons:
    merge_winter_nc(
        f"Runoff/Winter_2024_2025/runoff_dec_{dec_year}.nc",
        f"Runoff/Winter_2024_2025/runoff_jan_feb_{jf_year}.nc",
        f"Runoff/Winter_2024_2025/runoff_winter_{dec_year}_{jf_year}.nc",
    )

print("Done!")
