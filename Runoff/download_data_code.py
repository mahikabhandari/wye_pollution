"""Download EFAS runoff_water_equivalent data for seasons over a specified area and extract the NetCDF file from the downloaded zip archive."""

import os
import zipfile

import cdsapi

dataset = "efas-historical"
request = {
    "system_version": ["version_5_0"],
    "variable": ["runoff_water_equivalent"],
    "model_levels": "surface_level",
    "hyear": ["2024"],
    "hmonth": ["12"],
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

client = cdsapi.Client()
client.retrieve(dataset, request).download(
    "Runoff/Winter_2024_2025/runoff_dec_2024.zip"
)


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


zip_path = "Runoff/Winter_2024_2025/runoff_dec_2024.zip"
output_nc = "Runoff/Winter_2024_2025/runoff_dec_2024.nc"

extract_nc_from_zip(zip_path, output_nc)
