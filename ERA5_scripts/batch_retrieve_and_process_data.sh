#!/bin/sh

# retrieve the data for 1979
# soil temperatures
for year in {1979..2020}; do
    for level in 1 2 3 4; do 
        sbatch -o log/retrieve_${level}_${year}.log -e log/retrieve_${level}_${year}.err --open-mode=truncate \
               -t 23:59:59 --wrap="/usr/bin/time -v python retrieve_era5_hourly.py --area 58.5 123.0 52.0 134.0 stl$level $year" 
    done
done

# other instantaneous fields
for year in {1979..2020}; do
    for var in sd msl d2m t2m u10 v10; do 
        sbatch -o log/retrieve_${var}_${year}.log -e log/retrieve_${var}_${year}.err --open-mode=truncate \
               -t 23:59:59 --wrap="/usr/bin/time -v python retrieve_era5_hourly.py --area 58.5 123.0 52.0 134.0 $var $year"
    done
done

# hourly aggregates - min and max temperatures and precipitation
for year in {1979..2020}; do
    for var in mx2t mn2t tp; do 
        sbatch -o log/retrieve_${var}_${year}.log -e log/retrieve_${var}_${year}.err --open-mode=truncate \
               -t 23:59:59 --wrap="/usr/bin/time -v python retrieve_era5_hourly.py --area 58.5 123.0 52.0 134.0 $var $year" 
    done
done

# for aggregated hourly data, read the data for Jan 1980
# this is because the time coordinate in the CDS is the time at the *end* of the hour
# so requesting data for 1979 (that is data from 00:00 on 1 Jan 1979 to 23:00 on 31 Dec 1979)
# actually retrieves the final field for 1978 (well, it would if there were data for 1979)
# but omits the final field for 1979 (as its time coordinate is 00:00 on 1 Jan 1980, i.e. not in 1979)
for var in mx2t mn2t tp; do python retrieve_era5_hourly.py --area 58.5 123.0 52.0 134.0 $var 2021 01; done

# process the data to make hourly, 3-hourly and daily netCDF files
# note: for aggregated data (identified by entries in era5_variables.txt) this reads data both
#       from the specified year and the following year to make annual files that really do contain
#       all the data for the requested year
for year in {1979..2020}; do
    for dur in hourly 3-hourly daily; do
        for level in 1 2 3 4; do python process_era5_hourly_from_grib.py --area 64.0 127.5 57.5 141.5 --duration $dur stl${level} $year; done
        for var in mx2t mn2t tp; do python process_era5_hourly_from_grib.py --area 64.0 127.5 57.5 141.5 --duration $dur $var $year; done
        for var in sd msl d2m t2m u10 v10; do python process_era5_hourly_from_grib.py --area 64.0 127.5 57.5 141.5 --duration $dur $var $year; done
    done
done

