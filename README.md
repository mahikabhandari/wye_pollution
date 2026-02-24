# My MSci Project 
## Requirements

This requires the custom `funmixer` library available from https://github.com/AlexLipp/funmixer. Other dependencies are listed in the `environment.yaml` file. You can create a conda environment to run these scripts with all the other dependencies installed using:

```bash
conda env create -f environment.yaml
conda activate wye-inversion
```

The preprocessing script requires a D8 raster that must be downloaded separately from: https://zenodo.org/records/14238014.

## Contents 

- Preprocessing Python Script (`preprocessing_newd8.py`)
- Input data
   - Sites.csv
   - Phosphate_date.csv