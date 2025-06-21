# reanalysis_validation
Repository for downloading, processing and validating ERA5 reanalysis data against met station observations.

Use this code to download ERA5 data for a specified polygon from the CDS database, using the `cdsapi` and `iris` packages, and convert from grib to netCDF. Co-ordinates below the equator and/or west of the prime meridian are negative. The code can easily be adapted to download ERA-Interim or ERA5-Land data: just ensure the target variables and download path match those from the CDS database.

Order of processiong (files found in the `ERA5_scripts` folder):
* Ensure the full list of variables is correct in `era5_variables.txt` and define the target variables to be downloaded in `era5_grib_variable_info.txt`.
* Run the `retrieve_era5_hourly.py` script using the batch command located in `batch_retrieve_and_process_data.sh`.
* Process the data using `process_era5_hourly_from_grib.py` script using the batch command also located in the shell file.
* Additional processing options can be found in `process_era5_6-hourly_from_grib.py`.

Then extract the data at a single, specified point, related to the coordinates of the meteorological station using the `Extract_data_at_point.ipynb` Jupyter Notebook. Data will be interpolated to the exact point from the original coarse resolution. Download the meteorological station data separately and store in a csv file.

Finally, analyse the data either at station-level (`Analyse_data_at_single_station.ipynb`) or across all regions and stations (`Analyse_data_all_stations.ipynb`), and make plots.

Note that soil temperature data are treated separately as there are often large gaps in the soil temp data from met stations.

Retrieval and processing code correct as of Summer 2023.
