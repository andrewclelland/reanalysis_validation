#!/bin/sh

# retrieve the data for 1979
# soil temperatures
for level in 1 2 3 4; do python retrieve_era5-land_hourly.py --area 64.0 83.0 58.5 95.0 stl${level} 1979; done

# other instantaneous fields
for var in sd msl d2m t2m u10 v10; do python retrieve_era5-land_hourly.py --area 64.0 83.0 58.5 95.0 $var 1979; done

# hourly aggregates - total precipitation
for var in tp; do python retrieve_era5-land_hourly.py --area 64.0 83.0 58.5 95.0 $var 1979; done

# for aggregated hourly data, read the data for Jan 1980
# this is because the time coordinate in the CDS is the time at the *end* of the hour
# so requesting data for 1979 (that is data from 00:00 on 1 Jan 1979 to 23:00 on 31 Dec 1979)
# actually retrieves the final field for 1978 (well, it would if there were data for 1979)
# but omits the final field for 1979 (as its time coordinate is 00:00 on 1 Jan 1980, i.e. not in 1979)
for var in mx2t mn2t tp; do python retrieve_era5-land_hourly.py --area 64.0 83.0 58.5 95.0 $var 1980 01; done


# process the data to make hourly, 3-hourly and daily netCDF files
# note: for aggregated data (identified by entries in era5_variables.txt) this reads data both
#       from the specified year and the following year to make annual files that really do contain
#       all the data for the requested year
for dur in hourly 3-hourly daily; do
    for level in 1 2 3 4; do python process_era5-land_hourly_from_grib.py --area 64.0 83.0 58.5 95.0 --duration $dur stl${level} 1979; done
    for var in tp; do python process_era5-land_hourly_from_grib.py --area 64.0 83.0 58.5 95.0 --duration $dur $var 1979; done
    for var in sd msl d2m u10 v10; do python process_era5-land_hourly_from_grib.py --area 64.0 83.0 58.5 95.0 --duration $dur $var 1979; done
done

# handle t2m specially
# at the hourly timescale we only have t2m itself
python process_era5-land_hourly_from_grib.py --area 64.0 83.0 58.5 95.0 --duration hourly t2m 1979

# at 3-hourly and daily timescales we have mean, min and max of hourly t2m
for dur in 3-hourly daily; do
    for agg in mean minimum maximum; do
        python process_era5-land_hourly_from_grib.py --area 64.0 83.0 58.5 95.0 --duration $dur --aggregation $agg t2m 1979
    done
done
