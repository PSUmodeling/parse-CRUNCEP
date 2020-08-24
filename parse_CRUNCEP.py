#!/usr/bin/env python3

'''
Script to parse CRUNCEP data to generate Cycles weather files.
Usage: parse_CRUNCEP.py YYYY YYYY DATA_PATH"
'''

import math
import sys
import calendar
import numpy as np
from math import log
from statistics import mean
from datetime import timedelta, date, datetime
from netCDF4 import Dataset

def Closest(lat, lon, nc):

    best_y = (np.abs(nc.variables['lat'][:] - lat)).argmin()
    best_x = (np.abs(nc.variables['lon'][:] - lon)).argmin()

    return (best_y, best_x, nc['lat'][best_y], nc['lon'][best_x])


def read_var(y, x, k, tpqwl, solr, prec):

    _prcp  = prec['PRECTmms'][k, y, x]
    _temp  = tpqwl['TBOT'][k, y, x]
    _wind  = tpqwl['WIND'][k, y, x]
    _solar = solr['FSDS'][k, y, x]
    _pres  = tpqwl['PSRF'][k, y, x]
    _spfh  = tpqwl['QBOT'][k, y, x]

    es = 611.2 * math.exp(17.67 * (_temp - 273.15) / (_temp - 273.15 + 243.5))
    ws = 0.622 * es / (_pres - es)
    w = _spfh / (1.0 - _spfh)
    _rh = w / ws
    if _rh > 1.0:
        _rh = 1.0

    return (_prcp, _temp, _wind, _solar, _rh, _pres)

def satvp(temp):

    return .6108 * math.exp(17.27 * temp / (temp + 237.3))


def ea(patm, q):

    return patm * q / (0.622 * (1 - q) + q)


def process_month(year, month, y, x, path, data, pres):

    '''
    Process one month of CRUNCEP data and convert it to Cycles input
    '''

    print('%4.4d-%2.2d' % (year, month))

    first_day = 1

    # Temperature, pressure, humidity, wind
    fn = '%s/clmforc.cruncep.V7.c2016.0.5d.TPQWL.%4.4d-%2.2d.nc' % (path,
                                                                    year,
                                                                    month)
    tpqwl = Dataset(fn, 'r')

    # Solar radiation
    fn = '%s/clmforc.cruncep.V7.c2016.0.5d.Solr.%4.4d-%2.2d.nc' % (path,
                                                                   year,
                                                                   month)
    solr = Dataset(fn, 'r')

    # Precipitation
    fn = '%s/clmforc.cruncep.V7.c2016.0.5d.Prec.%4.4d-%2.2d.nc' % (path,
                                                                   year,
                                                                   month)
    prec = Dataset(fn, 'r')

    # Get number of days in current file. CRUNCEP data are 6-hourly
    ndays = int(len(tpqwl.variables['time']) / 4)

    # Number of locations of interest
    nloc = len(y)

    for d in range(ndays):
        # Current day
        t = datetime(year, month, 1) + timedelta(days = d)

        # Initialize daily variables
        prcp = [0.0] * nloc
        tx = [-999.0] * nloc
        tn = [999.0] * nloc
        wind = [0.0] * nloc
        solar = [0.0] * nloc
        rhx = [-999.0] * nloc
        rhn = [999.0] * nloc
        counter = 0

        for ktime in range(d * 3, d * 3 + 3):
            for kloc in range(nloc):
                (_prcp, _temp, _wind, _solar, _rh, _pres) = read_var(y[kloc],
                                                                     x[kloc],
                                                                     ktime,
                                                                     tpqwl,
                                                                     solr,
                                                                     prec)

                prcp[kloc] += _prcp

                if _temp > tx[kloc]:
                    tx[kloc] = _temp

                if _temp < tn[kloc]:
                    tn[kloc] = _temp

                wind[kloc] += _wind

                solar[kloc] += _solar

                if _rh > rhx[kloc]:
                    rhx[kloc] = _rh

                if _rh < rhn[kloc]:
                    rhn[kloc] = _rh

                pres[kloc].append(float(_pres))

            counter += 1

        for kloc in range(nloc):
            prcp[kloc] /= float(counter)
            prcp[kloc] *= 86400.0

            wind[kloc] /= float(counter)

            solar[kloc] /= float(counter)
            solar[kloc] *= 86400.0 / 1.0E6

            rhx[kloc] *= 100.0
            rhn[kloc] *= 100.0

            tx[kloc] -= 273.15
            tn[kloc] -= 273.15

            data[kloc] = data[kloc] + \
                '%-16s%-8.4f%-8.2f%-8.2f%-8.4f%-8.2f%-8.2f%-8.2f\n' \
                %(t.strftime('%Y    %j'), prcp[kloc], tx[kloc], tn[kloc],
                  solar[kloc], rhx[kloc], rhn[kloc], wind[kloc])

        first_day = 0

    # CRUNCEP does not have leap years. Copy Feb 28th data for Feb 29th
    if (calendar.isleap(year) and month == 2):
        t = datetime(year, 2, 29)
        for kloc in range(nloc):
            data[kloc] = data[kloc] + \
                '%-16s%-8.4f%-8.2f%-8.2f%-8.4f%-8.2f%-8.2f%-8.2f\n' \
                %(t.strftime('%Y    %j'), prcp[kloc], tx[kloc], tn[kloc],
                  solar[kloc], rhx[kloc], rhn[kloc], wind[kloc])

    tpqwl.close()
    prec.close()
    solr.close()

    return data


def main():

    if (len(sys.argv) != 4):
        print("Illegal number of parameters.")
        print("Usage: parse_CRUNCEP.py YYYY YYYY DATA_PATH")
        sys.exit(0)

    start_year = int(sys.argv[1])
    end_year = int(sys.argv[2])
    data_path = sys.argv[3]

    filepath = 'location.txt'
    cnt = 0
    y = []
    x = []
    outfp = []

    # Read location file
    with open(filepath) as fp:
        for _, line in enumerate(fp):
            li=line.strip()
            if not (li.startswith("#") or li.startswith("L")):
                cnt += 1

                nums = line.split()
                lat = float(nums[0])
                lon = float(nums[1])

                print('Processing data for {0}, {1}'.format(lat, lon))

                if lon < 0.0:
                    lon += 360.0

                # Get nearest grid
                nc_path = '%s/clmforc.cruncep.V7.c2016.0.5d.TPQWL.%4.4d-01.nc' \
                    %(data_path, start_year)

                nc = Dataset(nc_path, 'r')

                (_y, _x, grid_lat, grid_lon) = Closest(lat, lon, nc)
                x.append(_x)
                y.append(_y)

                nc.close()

                # Create output files, i.e., Cycles weather files
                if grid_lat < 0.0:
                    lat_str = '%.2fS' %(abs(grid_lat))
                else:
                    lat_str = '%.2fN' %(abs(grid_lat))

                if grid_lon > 180.0:
                    lon_str = '%.2fW' %(360.0 - grid_lon)
                else:
                    lon_str = '%.2fE' %(grid_lon)

                fname = 'cruncep' + lat_str + 'x' + lon_str + '.txt'
                outfp.append(open(fname, 'w'))
                outfp[cnt - 1].write('#lat %s, lon %s\n' % (nums[0], nums[1]))
                outfp[cnt - 1].write('LATITUDE            %s\n' % (nums[0]))

    data = []

    [data.append('YEAR    DOY     PP      TX      TN      SOLAR   RHX     RHN     WIND\n') for kloc in range(cnt)]

    pres = []
    [pres.append([]) for kloc in range(cnt)]

    # Loop through years and months
    for year in range(start_year, end_year + 1):
        for month in range(1, 12 + 1):
            process_month(year, month, y, x, data_path, data, pres)

    for kloc in range(cnt):
        outfp[kloc].write('ALTITUDE            %-.2f\n' \
            % (-8200.0 * log(mean(pres[kloc]) / 101325.0)))
        outfp[kloc].write('SCREENING_HEIGHT    2.0\n')

    [outfp[kloc].write(data[kloc]) for kloc in range(cnt)]

    [outfp[kloc].close() for kloc in range(cnt)]


main()

