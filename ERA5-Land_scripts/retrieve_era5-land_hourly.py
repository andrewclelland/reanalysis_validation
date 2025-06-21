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

def retrieve_era5_land_hourly(var_name, year, month, area, data_format):
    # get the CDS variable name and see if it is single level or pressure level data
    era5_variables = pd.read_table('era5_variables.txt', index_col='var_name')
    var_details = era5_variables.loc[var_name]
    cds_var = var_details.CDS_name
    cds_dataset = 'reanalysis-era5-land'

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

    month_str = '{:02d}'.format(month) if month is not None else ''
    target = "{}/era5-land_hourly_{}_{:04d}{}.{}".format(download_dir, var_name, year, month_str, extension)

    request = {'variable': cds_var,
               'year': '{:04d}'.format(year),
               'time': ['{:02d}:00'.format(hour) for hour in range(0, 24)],
               'day': ['{:02d}'.format(day) for day in range(1, 32)],
               'format': data_format}

    request['month'] = '{:02d}'.format(month) if month is not None \
                       else ['{:02d}'.format(month) for month in range(1, 13)]

    if area is not None:
        request['area'] = area

    print('Performing retrieval:\n\nCDS data set: {}\nRequest: {}\nTarget: {}\n'.format(cds_dataset, request, target))

    c.retrieve(cds_dataset, request, target)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--format", default="grib", help="Data format (grib/netcdf)")
    parser.add_argument("--area", default=None, nargs=4, help="Area (N W S E)", type=float)
    parser.add_argument("var", help="NetCDF variable name")
    parser.add_argument("year", help="year", type=int)
    parser.add_argument("month", help="Month", type=int, nargs='?', default=None)
    args = parser.parse_args()

    print('At {}, started'.format(time.strftime('%Y-%m-%d %H:%M:%S')))

    retrieve_era5_land_hourly(args.var, args.year, args.month, args.area, args.format)

    print('At {}, finished'.format(time.strftime('%Y-%m-%d %H:%M:%S')))
