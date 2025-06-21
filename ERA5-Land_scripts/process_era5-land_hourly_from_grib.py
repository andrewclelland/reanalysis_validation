import os
import pandas as pd
import argparse
import iris
from iris.time import PartialDateTime
from iris.coords import CellMethod
import iris.coord_categorisation as icc
import iris_grib
import iris_grib._grib_cf_map as grcf
import numpy as np
import scipy.constants
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import time


# add CF mappings for additional GRIB variables
grcf.GRIB1_LOCAL_TO_CF[grcf.G1LocalParam(1, 128, 98, 133)] = grcf.CFName('specific_humidity', None, 'kg kg-1')
grcf.GRIB1_LOCAL_TO_CF[grcf.G1LocalParam(1, 128, 98, 134)] = grcf.CFName('surface_air_pressure', None, 'Pa')
grcf.GRIB1_LOCAL_TO_CF[grcf.G1LocalParam(1, 128, 98, 137)] = grcf.CFName('atmosphere_mass_content_of_water_vapor', None, 'kg m-2')
grcf.GRIB1_LOCAL_TO_CF[grcf.G1LocalParam(1, 128, 98, 139)] = grcf.CFName(None, 'soil temperature level 1', 'K')
grcf.GRIB1_LOCAL_TO_CF[grcf.G1LocalParam(1, 128, 98, 144)] = grcf.CFName('lwe_thickness_of_snowfall_amount', None, 'm')
grcf.GRIB1_LOCAL_TO_CF[grcf.G1LocalParam(1, 128, 98, 170)] = grcf.CFName(None, 'soil temperature level 2', 'K')
grcf.GRIB1_LOCAL_TO_CF[grcf.G1LocalParam(1, 128, 98, 172)] = grcf.CFName('land_binary_mask', None, '1')
grcf.GRIB1_LOCAL_TO_CF[grcf.G1LocalParam(1, 128, 98, 182)] = grcf.CFName('lwe_thickness_of_water_evaporation_amount', None, 'm')
grcf.GRIB1_LOCAL_TO_CF[grcf.G1LocalParam(1, 128, 98, 183)] = grcf.CFName(None, 'soil temperature level 3', 'K')
grcf.GRIB1_LOCAL_TO_CF[grcf.G1LocalParam(1, 128, 98, 201)] = grcf.CFName(None, 'maximum temperature at 2 metres since previous post-processing', 'K')
grcf.GRIB1_LOCAL_TO_CF[grcf.G1LocalParam(1, 128, 98, 202)] = grcf.CFName(None, 'minimum temperature at 2 metres since previous post-processing', 'K')
grcf.GRIB1_LOCAL_TO_CF[grcf.G1LocalParam(1, 128, 98, 228)] = grcf.CFName('lwe_thickness_of_precipitation_amount', None, 'm')
grcf.GRIB1_LOCAL_TO_CF[grcf.G1LocalParam(1, 128, 98, 236)] = grcf.CFName(None, 'soil temperature level 4', 'K')

iris_grib.grib_phenom_translation._GRIB1_CF_TABLE = iris_grib.grib_phenom_translation._make_grib1_cf_table()


DEFAULT_CHUNK_SIZE = 4*1024*1024

def default_chunk_sizes(cube):
    sizes = {}
    for coord in cube.coords(dim_coords=True):
        # code derived from lines 348:394 of https://github.com/Unidata/netcdf-c/blob/v4.6.1/libsrc4/nc4var.c
        suggested_size = int(((DEFAULT_CHUNK_SIZE / cube.data.nbytes) **
                          (1.0 / len(cube.shape))) * len(coord.points) - 0.5)
        if suggested_size < 1:
            suggested_size = 1
        elif suggested_size > len(coord.points):
            suggested_size = len(coord.points)
        num_chunks = (len(coord.points) + suggested_size - 1) // suggested_size        
        if num_chunks > 0:
            overhang = (num_chunks * suggested_size) - len(coord.points)
            suggested_size -= overhang // num_chunks
        sizes[coord.name()] = suggested_size

    return sizes

def format_longitude(lon):
    return '{}{}'.format(abs(lon), '' if lon == 0 else 'E' if lon > 0 else 'W')

def format_latitude(lat):
    return '{}{}'.format(abs(lat), '' if lat == 0 else 'N' if lat > 0 else 'S')

def process_era5_land_hourly_from_grib(var_name, area, duration, aggregation,
                                       start_year, end_year, start_month, end_month):
    # variable types and aggregation methods
    time_at_end_var_types = ['accumulated', 'mean rate', 'minimum', 'maximum']
    aggregations = {'accumulated': 'total', 'instantaneous': 'mean',
                    'invariant': 'mean', 'mean rate': 'mean',
                    'minimum': 'minimum', 'maximum': 'maximum'}
    aggregators = {'total': iris.analysis.SUM, 'mean': iris.analysis.MEAN,
                    'minimum': iris.analysis.MIN, 'maximum': iris.analysis.MAX}
    agg_prefixes = {'total': 'total', 'mean': 'mean', 'minimum': 'min', 'maximum': 'max'}

    # find out whether the variable is instantaneous or accumulated/rate/min/max
    era5_variables = pd.read_table('era5_variables.txt', index_col='var_name')
    var_details = era5_variables.loc[var_name]
    var_type = var_details.type
    time_at_end = var_type in time_at_end_var_types

    # see whether the aggregation is specified and use the default if not
    default_aggregation = (aggregation == 'default')
    if default_aggregation:
        aggregation = aggregations[var_type]

    # map the variable to long name, standard name and cell method
    era5_grib_variable_info = pd.read_table('era5_grib_variable_info.txt', index_col='var_name',
                                            keep_default_na=False)
    try:
        var_details = era5_grib_variable_info.loc[var_name]
        long_name = var_details.long_name
        standard_name = var_details.standard_name
        cell_method = var_details.cell_method
    except KeyError:
        long_name, standard_name, cell_method = None, None, None

    # data directory
    data_dir = []
    if area is not None:
        data_dir.append('{s}-{n}_{w}-{e}'.format(
                        n=format_latitude(area[0]),
                        w=format_longitude(area[1]),
                        s=format_latitude(area[2]),
                        e=format_longitude(area[3])))
    data_dir = '/'.join(data_dir)

    # get the set of yearly or monthly filenames spanning the requested date range
    download_dir = '../{}/download'.format(data_dir)
    fn_pattern = '{}/era5-land_hourly_{}_{:04d}*.grib'
    filenames = []
    for year in range(start_year, end_year+1):
        filenames.append(fn_pattern.format(download_dir, var_name, year))

    # if we are reading accumulated/rate/min/max data and the date range ends at the end of a year
    # add the filenames for the month or year after the final year
    # (as this will have the data for the final time step)
    if end_month == 12 and time_at_end:
        final_month = datetime(end_year, end_month, 1) + relativedelta(months=1)
        filenames.append(fn_pattern.format(download_dir, var_name, final_month.year))

    # read the hourly data
    print('Reading data from {}'.format(filenames))
    hourly_cube = iris.load(filenames).concatenate_cube()

    # change the time origin to be 00Z on 01-Jan-1900
    time_coord = hourly_cube.coord(axis='t')
    time_coord.convert_units(time_coord.units.origin.replace('1970', '1900'))

    # identify the spatial coordinate system (spherical Earth, radius 6367.47 km)
    # see https://confluence.ecmwf.int/display/CKB/ERA5%3A+What+is+the+spatial+reference
    cs = iris.coord_systems.GeogCS(6367470.0)
    hourly_cube.coord(axis='x').coord_system = cs
    hourly_cube.coord(axis='y').coord_system = cs

    # extract the specified range of dates
    time_coord_name = time_coord.name()
    start_date_pdt = PartialDateTime(start_year, start_month)
    start_con = iris.Constraint(coord_values=
                                {time_coord_name:lambda cell: start_date_pdt <= cell.point})
    end_date_pdt = PartialDateTime(end_year, end_month)
    end_con = iris.Constraint(coord_values=
                              {time_coord_name:lambda cell: end_date_pdt >= cell.point})
    hourly_cube = hourly_cube.extract(start_con & end_con)

    # set the cube name details
    hourly_cube.var_name = var_name
    if long_name is not None:
        hourly_cube.long_name = long_name
    if standard_name is not None:
        hourly_cube.standard_name = standard_name
    if cell_method is not None and len(cell_method) != 0:
        hourly_cube.add_cell_method(iris.coords.CellMethod(cell_method, time_coord_name, '1 hour'))

    # aggregate the data if required
    long_name = hourly_cube.long_name
    if duration == 'hourly':
        output_cube = hourly_cube
    elif duration in ['3-hourly', 'daily']:
        # add categorical coordinates so the data can be aggregated and
        # identify the requierd duration (so partial aggregates can be removed)
        time_coord = hourly_cube.coord(axis='t')
        icc.add_year(hourly_cube, time_coord)
        icc.add_month_number(hourly_cube, time_coord, 'month')
        icc.add_day_of_month(hourly_cube, time_coord)
        agg_coords = ['year', 'month', 'day_of_month']
        required_duration = 24
        if duration == '3-hourly':
            icc.add_categorised_coord(hourly_cube, '3-hourly_period', time_coord,
                                      lambda coord, value: time_coord.units.num2date(value).hour // 3)
            agg_coords.append('3-hourly_period')
            required_duration = 3

        # aggregate the data
        output_cube = hourly_cube.aggregated_by(agg_coords, aggregators[aggregation])

        # rename the variable if the default aggregation was not used
        if not default_aggregation:
            output_cube.var_name = '{}_{}'.format(agg_prefixes[aggregation], var_name)
            output_cube.long_name = '{} {}'.format(aggregation, long_name)

        # for instantaneous data, the actual duration is one hour less than the specified
        # aggregation: for example, a mean from 00:00 to 23:00 has a duration of 23 hours
        if not time_at_end:
            required_duration -= 1

        # keep only aggregated fields with the required duration (to remove partial aggregates)
        tdelta = timedelta(hours=required_duration)
        time_coord_name = output_cube.coord(axis='t').name()
        dur_con = iris.Constraint(coord_values=
                                  {time_coord_name:lambda t: (t.bound[1] - t.bound[0]) == tdelta})
        output_cube = output_cube.extract(dur_con)
    else:
        raise ValueError('Unknown duration: {}'.format(duration))

    # cast the data to float
    output_cube.data = output_cube.data.astype(np.float32)

    # if the cube is masked, make all the masked values equal its fill value
    fill_value = None
    if np.ma.is_masked(output_cube.data):
        fill_value = output_cube.data.fill_value
        output_cube.data.data[output_cube.data.mask] = fill_value

    # add a title attribute to identify the data
    aggregation_desc = '' if duration == 'hourly' else '{} {} of '.format(duration, aggregation)
    output_cube.attributes['title'] = 'ERA5-Land reanalysis {}hourly {}'.format(aggregation_desc, long_name)

    # work out the filename for the cube
    odir = ['/work/scratch-pw/{}/era5'.format(os.getenv('USER'))]
    if data_dir != '':
        odir.append(data_dir)
    odir.extend([duration, output_cube.var_name])
    odir = '/'.join(odir)
    nc_start = 'era5-land_{}'.format(duration)
    if start_year != end_year:
        if start_month == 1 and end_month == 12:
            nc_base = '{}_{}_{:04d}-{:04d}'.format(nc_start,
                      output_cube.var_name, start_year, end_year)
        else:
            nc_base = '{}_{}_{:04d}{:02d}-{:04d}{:02d}'.format(nc_start,
                      output_cube.var_name, start_year, start_month, end_year, end_month)
    elif start_month != end_month:
        if start_month == 1 and end_month == 12:
            nc_base = '{}_{}_{:04d}'.format(nc_start, output_cube.var_name, start_year)
        else:
            nc_base = '{}_{}_{:04d}{:02d}-{:04d}{:02d}'.format(nc_start,
                      output_cube.var_name, start_year, start_month, end_year, end_month)
    else:
        nc_base = '{}_{}_{:04d}{:02d}'.format(nc_start,
                  output_cube.var_name, start_year, start_month)

    # create the directory if required
    if not os.path.exists(odir):
        os.makedirs(odir)

    # save the cube to a temporary (uncompressed) file
    tmp_nc_name = '{}/{}_tmp.nc'.format(odir, nc_base)
    print('Saving data to temporary file {}'.format(tmp_nc_name))
    iris.save(output_cube, tmp_nc_name)

    # compress the file
    nc_name = '{}/{}.nc'.format(odir, nc_base)
    chunkspec = ','.join(['{}/{}'.format(dim, chunk_size) for dim, chunk_size in 
                          default_chunk_sizes(output_cube).items()])
    compress_cmd = 'nccopy -d 9 -c {} {} {}'.format(chunkspec, tmp_nc_name, nc_name)
    print('Compressing file: {}'.format(compress_cmd))
    os.system(compress_cmd)
    os.remove(tmp_nc_name)

    # copy the output file from the work area to the appropriate place
    fdir = ['..']
    if data_dir != '':
        fdir.append(data_dir)
    fdir.extend([duration, output_cube.var_name])
    fdir = '/'.join(fdir)
    if not os.path.exists(fdir):
        os.makedirs(fdir)
    copy_cmd = 'rsync -tvz {} {}'.format(nc_name, fdir)
    print('Copying file by running: {}'.format(copy_cmd))
    os.system(copy_cmd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("var", help="NetCDF variable name")
    parser.add_argument("--area", default=None, nargs=4, help="Area (N W S E)", type=float)
    parser.add_argument("--duration", choices=["hourly", "3-hourly", "daily"], default="hourly")
    parser.add_argument("--aggregation", choices=["mean", "total", "minimum", "maximum"], default="default")
    parser.add_argument("year", help="Start year", type=int)
    parser.add_argument("end_year_or_start_month", help="End year or Start month", type=int, nargs='?', default=None)
    parser.add_argument("end_month", help="Month or End month", type=int, nargs='?', default=None)
    args = parser.parse_args()

    if args.end_year_or_start_month is None:
        start_year = args.year
        end_year = start_year
        start_month = 1
        end_month = 12
    elif args.end_year_or_start_month >= args.year:
        start_year = args.year
        end_year = args.end_year_or_start_month
        if args.end_month is None:
            start_month = 1
            end_month = 12
        else:
            start_month = args.end_month
            end_month = start_month
    else:
        start_year = args.year
        end_year = start_year
        if args.end_month is None:
            start_month = args.end_year_or_start_month
            end_month = args.end_year_or_start_month
        else:
            start_month = args.end_year_or_start_month
            end_month = args.end_month

    print('At {}, started'.format(time.strftime('%Y-%m-%d %H:%M:%S')))

    process_era5_land_hourly_from_grib(args.var, args.area, args.duration, args.aggregation,
                                       start_year, end_year, start_month, end_month)

    print('At {}, finished'.format(time.strftime('%Y-%m-%d %H:%M:%S')))
