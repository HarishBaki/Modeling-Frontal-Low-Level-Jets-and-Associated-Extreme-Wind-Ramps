# Modelling Frontal Low-Level Jets and Associated Extreme Wind Power Ramps over the North Sea

This repository contains scripts and data for analyzing wind turbine data using CERRA and ERA5 datasets. The necessary data is stored separately on Zenodo.
This repository serves as a data analysis for the article titled "Modeling Frontal Low-Level Jets and Associated Extreme Wind Power Ramps over the North Sea", published in the journal of Wind Energy Science. 

## Repository Structure

```
Analysis.ipynb
borssele_wind_turbines.nc
Download_CERRA_convert_ERA5.py
download_ERA5.py
myoutfields.txt
tslist
Vtable.CERRA
wind-turbine-1.tbl
wind-turbine-2.tbl
wind-turbine-3.tbl
wind-turbine-4.tbl
wind-turbine-5.tbl
windturbines.txt
data_processing/
    CERRA_and_ERA5_wind_power_through_chebyshev_coefifcients.ipynb
    compute_turbine_power.py
    extract_POWER.py
    libraries.py
namelist.input_files/
    ...
namelist.wps_files/
    ...
WES_dataset/
    CERRA/
    ERA5/
    Elia_dataset/
    Lidar_observations/
    WRF_simulations/
```

## Data

The necessary data for this analysis is stored separately on Zenodo. Please download the data at [zenodo](10.5281/zenodo.15033463) and place it in the appropriate directories.

## Namelist Files

The `namelist.wps_files` directory contains the namelist files required for WPS, and the `namelist.input_files` directory contains the namelist files required for WRF.

## Turbine Data

The locations of turbines, turbine power curves, and turbine specifications are provided through `.tbl` and `.txt` files:

- `wind-turbine-1.tbl`
- `wind-turbine-2.tbl`
- `wind-turbine-3.tbl`
- `wind-turbine-4.tbl`
- `wind-turbine-5.tbl`
- `windturbines.txt`

## Scripts

Detailed scripts are provided to download CERRA and ERA5 initial and boundary conditions:

- `Download_CERRA_convert_ERA5.py`
- `download_ERA5.py`

## Analysis

The entire analysis is done in the `Analysis.ipynb` notebook. This notebook contains all the steps and code required to perform the analysis.

## How to Run

1. Download the necessary data from Zenodo and place it in the appropriate directories.
2. Run the scripts to download CERRA and ERA5 initial and boundary conditions.
3. Open the `Analysis.ipynb` notebook and follow the steps to perform the analysis.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- The data used in this analysis is provided by Zenodo.
- The scripts and analysis are based on the work of the contributors to this repository.