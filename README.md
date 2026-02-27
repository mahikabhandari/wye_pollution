# My MSci Project 
## Requirements

This requires the custom `funmixer` library available from https://github.com/AlexLipp/funmixer. Other dependencies are listed in the `environment.yaml` file. You can create a conda environment to run these scripts with all the other dependencies installed using:

```bash
conda env create -f environment.yaml
conda activate wye-inversion
```

The preprocessing script requires a D8 raster that must be downloaded separately from: https://zenodo.org/records/14238014.

## Key Contents 

- Input data - `Data_PP/Original/`
- Preprocessing scripts for Wye citizen science phosphorous data - `Preprocessing_python_scripts/`
- EFAS Runoff data download, processing and outputs - `Runoff/`
- Scripts and outputs for phosphorous source inversion: `Modelling/` 

## Instructions

1. Run the preprocessing Python Script (`preprocessing_newd8.py`) -- This should be run first, as it prepares the input data for the subsequent analysis. 
2. Next, download the runoff data from EFAS by running `Runoff/download_data_code.py` -- Automatically downloads data from EFAS broken down by season
3. Reproject these data onto same grid as our D8, and average the seasonal runoff using `Runoff/runoff_average_reprojection.py`.
4. We next run `Modelling/Seasonal_phosphate_averages.py` to calculate the seasonal average Phosphate concentration for each subcatchment then calculate the average runoff in each catchment using `Runoff/runoff_catchment.py`.
5. Now we are ready to invert the phosphorous data using `Modelling/automated_seasons.py`