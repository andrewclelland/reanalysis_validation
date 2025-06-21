# reanalysis_validation
Repository for downloading, processing and validating ECMWF reanalysis data against met station observations.

Use this code to download ERA-Interim/ERA5/ERA5-Land (change as appropriate) data for a specified polygon from the CDS database, using the `cdsapi` package, and convert from grib to netCDF. North and East co-ordinates are positive; South and West are negative.

Then extract the data at a single, specified point, related to the coordinates of the meteorological station, using the `iris` package. Data will be interpolated to the exact point from the original coarse resolution. Download the meteorological station data separately and store in a csv file.

Finally, analyse the data either at station-level or across all regions and stations, and make plots.
