"""
Code to retrieve the data from the CDS server. Ensure the chosen variables are present in the 'era5_grib_variable_info.txt' file beforehand, and they match those from the whole list in 'era5_variables.txt'.

Also ensure the `cdsapi` package is installed and setup locally beforehand.
"""

import os
import argparse
import cdsapi
import pandas as pd
import time

c = cdsapi.Client()

def format_longitude(lon):
    return '{}{}'.format(abs(lon), '' if lon == 0 else 'E' if lon > 0 else 'W')

def format_latitude(lat):
    return '{}{}'.format(abs(lat), '' if lat == 0 else 'N' if lat > 0 else 'S')

def retrieve_era5_hourly(var_name, year, month, level, area, data_format):
    # CDS data sets
    cds_datasets = {'sfc': 'reanalysis-era5-single-levels',
                    'pl': 'reanalysis-era5-pressure-levels'}

    # get the CDS variable name and see if it is single level or pressure level data
    era5_variables = pd.read_table('era5_variables.txt', index_col='var_name')
    var_details = era5_variables.loc[var_name]
    cds_var = var_details.CDS_name
    levels = var_details.levels
    dataset = 'pl' if levels == 'pressure' else 'sfc'
    cds_dataset = cds_datasets[dataset]

    # download directory
    download_dir = ['..']
    if area is not None:
        download_dir.append('{s}-{n}_{w}-{e}'.format(
                            n=format_latitude(area[0]),
                            w=format_longitude(area[1]),
                            s=format_latitude(area[2]),
                            e=format_longitude(area[3])))
    download_dir.append('download')
    download_dir = '/'.join(download_dir)

    # create the directory if required
    if not os.path.exists(download_dir):
        os.makedirs(download_dir)

    # filename extension
    extensions = {'netcdf':'nc', 'grib':'grib'}
    extension = extensions[data_format]

    var_with_level = var_name if dataset == 'sfc' else '{}_{}hPa'.format(var_name, level)
    month_str = '{:02d}'.format(month) if month is not None else ''
    target = "{}/era5_hourly_{}_{:04d}{}.{}".format(download_dir, var_with_level, year, month_str, extension)

    request = {'product_type': 'reanalysis',
               'variable': cds_var,
               'year': '{:04d}'.format(year),
               'time': ['{:02d}:00'.format(hour) for hour in range(0, 24)],
               'day': ['{:02d}'.format(day) for day in range(1, 32)],
               'format': data_format}

    request['month'] = '{:02d}'.format(month) if month is not None \
                       else ['{:02d}'.format(month) for month in range(1, 13)]

    if dataset == 'pl':
        request['pressure_level'] = [level]

    if area is not None:
        request['area'] = area

    print('Performing retrieval:\n\nCDS data set: {}\nRequest: {}\nTarget: {}\n'.format(cds_dataset, request, target))

    c.retrieve(cds_dataset, request, target)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--format", default="grib", help="Data format (grib/netcdf)")
    parser.add_argument("--level", default=None, help="Level (in hPa) - for pressure level data only", type=int)
    parser.add_argument("--area", default=None, nargs=4, help="Area (N W S E)", type=float)
    parser.add_argument("var", help="NetCDF variable name")
    parser.add_argument("year", help="year", type=int)
    parser.add_argument("month", help="Month", type=int, nargs='?', default=None)
    args = parser.parse_args()

    print('At {}, started'.format(time.strftime('%Y-%m-%d %H:%M:%S')))

    retrieve_era5_hourly(args.var, args.year, args.month, args.level, args.area, args.format)

    print('At {}, finished'.format(time.strftime('%Y-%m-%d %H:%M:%S')))
