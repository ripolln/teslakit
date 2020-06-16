#!/usr/bin/env python
# -*- coding: utf-8 -*-

#common
import time
import sys
import os
import os.path as op
from datetime import datetime, timedelta
import threddsclient
import signal

# pip
import numpy as np
import netCDF4 as nc4
import xarray as xr
from geographiclib.geodesic import Geodesic

# TODO:
#    - try/catch archivos incompletos
#    - multithreading

def geo_distance(lat1, lon1, lat2, lon2):
    'Returns geodesic distance between points in degrees'

    return Geodesic.WGS84.Inverse(
        lat1, lon1, lat2, lon2
    )['a12']

def join_ncs(p_ncs_folder, p_out_nc):
    '''
    Join .nc files downloaded from CSIRO
    p_ncs_folder - folder containing CSIRO downloaded netCDF4 files
    p_out_nc     - path to output netcdf file
    '''

    chk = 12 # n files for packs 

    # clean previous packs
    l_packs_del = sorted(
        [op.join(p_ncs_folder,n) for n in os.listdir(p_ncs_folder) if
         n.startswith('xds_packed_')])
    for f in l_packs_del:
        os.remove(f)

    # get files to join
    l_ncs = sorted(
        [op.join(p_ncs_folder,n) for n in os.listdir(p_ncs_folder) if
         n.endswith('.nc')])
    print('joining {0} netCDF4 files...'.format(len(l_ncs)), end='')

    # first join by chk size in .nc packs 
    c = 1
    while l_ncs:
        pack = l_ncs[:chk]
        l_ncs = l_ncs[chk:]

        xds_pack = xr.open_mfdataset(pack)
        p_pack = op.join(p_ncs_folder, 'xds_packed_{0:03d}.nc'.format(c))
        xds_pack.to_netcdf(p_pack,'w')
        xds_pack.close()
        c+=1

    # join pakcs
    l_packs = sorted(
        [op.join(p_ncs_folder,n) for n in os.listdir(p_ncs_folder) if
         n.startswith('xds_packed_')])
    xds_out = xr.open_mfdataset(l_packs)

    # uncompressed
    xds_out.to_netcdf(p_out_nc,'w')

    # compressed 
    #comp = dict(zlib=True, complevel=5)
    #encoding = {var: comp for var in xds_out.data_vars}
    #xds_out.to_netcdf(p_out_nc, 'w', encoding=encoding)

    # clean packs
    for p in l_packs:
        os.remove(p)

    print(' done.')
    return xr.open_dataset(p_out_nc)

def signal_handler(signum, frame):
    raise Exception("Timed out!")

def generate_urls(switch_db='gridded', grid_code='glob_24m', time_limits=None):
    '''
    Generate URL list for downloading csiro gridded/spec data
    switch_db - gridded / spec
    grid_code - 'pac_4m', 'pac_10m', 'aus_4m', 'aus_10m', 'glob_24m'
    time_limits - optional tuple (time_start, time_end) with 'yyyymm' format

    returns list of urls with monthly CSIRO data
    '''

    # parameters
    url_catg = 'http://data-cbr.csiro.au/thredds/catalog/catch_all/CMAR_CAWCR-Wave_archive/'
    url_dodsC = 'http://data-cbr.csiro.au/thredds/dodsC/catch_all/CMAR_CAWCR-Wave_archive/'
    url_a = 'CAWCR_Wave_Hindcast_aggregate/'

    # mount URLS
    url_xml = '{0}{1}{2}/catalog.xml'.format(url_catg, url_a, switch_db)
    cat = threddsclient.read_url(url_xml)
    down_ncs = sorted([f.name for f in cat.flat_datasets()])

    # filter grid
    if switch_db == 'gridded':
        down_ncs = [x for x in down_ncs if grid_code in x]

    # filter times (optional)
    if time_limits != None:
        ix1 = [time_limits[0] in u for u in down_ncs].index(True)
        ix2 = [time_limits[1] in u for u in down_ncs].index(True)
        down_ncs = down_ncs[ix1:ix2+1]

    l_urls = ['{0}{1}{2}/{3}'.format(url_dodsC, url_a, switch_db, nn) for nn in down_ncs]

    # mount URLs
    #if switch_db == 'gridded':
    #    # get data available inside gridded catalog
    #    url_xml = '{0}{1}gridded/catalog.xml'.format(url_catg, url_a)
    #    cat = threddsclient.read_url(url_xml)
    #    down_ncs = sorted([f.name for f in cat.flat_datasets()])
    #    down_ncs = [x for x in down_ncs if grid_code in x]  # filter grid

    #    l_urls = ['{0}{1}gridded/{2}'.format(
    #        url_dodsC, url_a, nn) for nn in down_ncs]

    #elif switch_db == 'spec':
    #    # get data available inside spec catalog
    #    url_xml = '{0}{1}spec/catalog.xml'.format(url_catg, url_a)
    #    cat = threddsclient.read_url(url_xml)
    #    down_ncs = sorted([f.name for f in cat.flat_datasets()])

    #    l_urls = ['{0}{1}spec/{2}'.format(
    #        url_dodsC, url_a, nn) for nn in down_ncs]

    return l_urls

# SPECTRAL
def url_extract_spec(url, stations_info, p_save):
    '''
    Extract spectral stations data.

    url           - url to remote csiro spec monthly file
    p_save        - path to store extracted stations netCDF4 file
    stations_info - list of tuples with stations information
    (index, id, lon, lat)
    '''

    # TODO: chunk data (maybe needed if many stations)

    # get stations info
    stations_ix = [t[0] for t in stations_info]
    stations_ID = [t[1] for t in stations_info]

    # extract file from url
    with xr.open_dataset(url) as xds_u:

        # iterate stations 
        for s_ix, s_ID in zip(stations_ix, stations_ID):
            p_station = op.join(p_save, '{0}'.format(s_ID))
            p_nc = op.join(p_station, op.basename(url))

            if not op.isdir(p_station): os.makedirs(p_station)
            if op.isfile(p_nc): continue

            # extract station
            xds_e = xds_u.isel(station=s_ix).squeeze()
            xds_e.to_netcdf(p_nc, 'w')

def download_spec_stations(l_urls, p_dl, stations_info):
    '''
    Download CSIRO spectral stations (using id).
    Also genrates a log file

    l_urls - list of urls to remote CSIRO spec files
    p_dl   - path to store downloaded stations

    stations_info: list of tuples with stations information:
    (index, id, lon, lat)
    '''

    # get stations info
    stations_ID = [t[1] for t in stations_info]

    # ini log 
    log = ['CSIRO SPEC STATIONS DOWNLOAD LOG\n\n']
    for si in stations_info:
        log.append('station ID: {0}. longitude: {1:.2f}, latitude: {2:.2f}\n'.format(
            si[1], si[2], si[3])
        )
    print(''.join(log))
    log.append('\ndownload start: {0}\n'.format(str(datetime.now())))

    # download stations
    print('downloading CSIRO spec data... {0} files'.format(len(l_urls)))
    for u in l_urls:
        print(op.basename(u), ' ... ', end='')
        sys.stdout.flush()
        url_extract_spec(u, stations_info, p_dl)
        print('downloaded.')

    # join .nc files
    for s_ID in stations_ID:
        p_st_dl = op.join(p_dl, '{0}'.format(s_ID))
        p_ncf = op.join(p_dl, '..', 'station_{0}.nc'.format(s_ID))
        join_ncs(p_st_dl, p_ncf)

    # end log
    log.append('download end:   {0}\n'.format(str(datetime.now())))
    with open(op.join(p_dl, '..', 'info_spec.txt'), 'w') as fW:
        fW.writelines(log)

def download_spec(p_store, point_list, time_limits=None):
    '''
    Locate and download CSIRO spectral stations data, then stores it on netcdf format

    point_list - query, list of tuples (lon_point, lat_point)
    p_store    - path to store downloaded stations
    time_limits - optional tuple (time_start, time_end) with 'yyyymm' format
    '''

    # paths
    p_dl = op.join(p_store, 'spec.download')  # temp storage of monthly files
    if not op.isdir(p_dl): os.makedirs(p_dl)

    # Generate URL list 
    l_urls = generate_urls('spec', time_limits=time_limits)

    # get stations lat,lon from first file
    with nc4.Dataset(l_urls[0], 'r') as ff:
        slons = ff['longitude'][0,:]
        slats = ff['latitude'][0,:]
        sids = ff['station'][:]

    # find nearest stations ids
    stations_ix = []
    for lon_p, lat_p in point_list:
        stations_ix.append((np.sqrt(np.square(slats-lat_p)+np.square(slons-lon_p))).argmin())

        # TODO calculate geodesic distance (degree)
        #geo_dist = []
        #for lon_ps, lat_ps in zip(slons, slats):
        #     geo_dist.append(geo_distance(lat_ps, lon_ps, p_lat, p_lon))
        #geo_dist = np.asarray(geo_dist)

    # join some station information in a list of tuples
    stations_info = [(x, sids[x], slons[x], slats[x]) for x in stations_ix]

    # download selected stations
    download_spec_stations(l_urls, p_dl, stations_info)

def download_spec_area(p_store, lon1, lat1, lon2, lat2, time_limits=None):
    '''
    Locate and download CSIRO spectral stations data, then stores it on netcdf format

    stations are searched inside area:
    lon1, lat1 - bottom left
    lon2, lat2 - top right
    p_store    - path to store downloaded stations
    time_limits - optional tuple (time_start, time_end) with 'yyyymm' format
    '''

    # paths
    p_dl = op.join(p_store, 'spec.download')  # temp storage of monthly files
    if not op.isdir(p_dl): os.makedirs(p_dl)

    # Generate URL list 
    l_urls = generate_urls('spec', time_limits=time_limits)

    # find stations inside area 
    with nc4.Dataset(l_urls[0], 'r') as ff:
        slons = ff['longitude'][0,:]
        slats = ff['latitude'][0,:]
        sids = ff['station'][:]

        stations_ix = np.where(
            (slons >= lon1) & (slons <= lon2) & \
            (slats >= lat1) & (slats <= lat2)
        )[0]

    # join some station information in a list of tuples
    stations_info = [(x, sids[x], slons[x], slats[x]) for x in stations_ix]

    # download selected stations
    download_spec_stations(l_urls, p_dl, stations_info)


# GRIDDED
def url_extract_gridded_points(url, points_info, p_save):
    '''
    Extract spectral points data.

    url           - url to remote csiro gridded monthly file
    p_save        - path to store extracted points netCDF4 file
    points_info   - list of tuples with points information:
    (index, id, lon, lat)
    '''

    # TODO: this may get stuck. add try/catch/timed signal 

    # get points info
    points_ID = [t[1] for t in points_info]
    points_lon = [t[2] for t in points_info]
    points_lat = [t[3] for t in points_info]

    # extract file from url
    with xr.open_dataset(url) as xds_u:

        # iterate points 
        for p_ID, p_lon, p_lat in zip(points_ID, points_lon, points_lat):
            p_point = op.join(p_save, '{0}'.format(p_ID))
            p_nc = op.join(p_point, op.basename(url))

            if not op.isdir(p_point): os.makedirs(p_point)
            if op.isfile(p_nc): continue

            # extract point
            xds_e = xds_u.sel(longitude=p_lon, latitude=p_lat).squeeze()
            xds_e.to_netcdf(p_nc, 'w')

def url_extract_gridded_area(url, ix1, iy1, ix2, iy2, p_save):
    '''
    Extract spectral points data.

    url             - url to remote csiro gridded monthly file
    p_save          - path to store extracted points netCDF4 file
    ix1,iy1,ix2,iy2 - nc file lon lat ini end indexes
    '''

    p_nc = op.join(p_save, op.basename(url))
    if op.isfile(p_nc): return

    # TODO: this may get stuck. add try/catch/timed signal 

    # extract file from url
    with xr.open_dataset(url) as xds_u:

        # extract area
        xds_e = xds_u.isel(
            longitude = slice(ix1, ix2+1),
            latitude = slice(iy1, iy2+1)
        )
        xds_e.to_netcdf(p_nc, 'w')

    #def getrfile(u, p_nc):
    #    'download unstable gridded files'
    #    min_timeout = 30

    #    while not op.isfile(p_nc):
    #        try:
    #            # we limit time available for operation
    #            signal.signal(signal.SIGALRM, signal_handler)
    #            signal.alarm(min_timeout*60)   # min*60 seconds

    #            # read file from url
    #            with xr.open_dataset(u) as xds_u:
    #                xds_temp = xds_u.isel(
    #                    longitude=slice(idx1,idx2+1),
    #                    latitude=slice(idy1,idy2+1))
    #                # save temp file
    #                xds_temp.to_netcdf(p_nc,'w')

    #        except Exception:
    #            # clean failed download and retry
    #            if op.isfile(p_nc):
    #                os.remove(p_nc)
    #            print('timed out. retry... ',end=' ')
    #            sys.stdout.flush()
    #        except:
    #            # clean failed download and retry
    #            if op.isfile(p_nc):
    #                os.remove(p_nc)
    #            print('failed. retry... ',end=' ')
    #            sys.stdout.flush()

    #        finally:
    #            signal.alarm(0)

def download_gridded_points(l_urls, p_dl, points_info):
    '''
    Download CSIRO gridded points to individual files.
    Also generates a log file

    l_urls - list of urls to remote CSIRO spec files
    p_dl   - path to store downloaded stations

    points_info: list of tuples with points information:
    (_, name, lon, lat)
    '''

    # get points info
    points_ID = [t[1] for t in points_info]

    # ini log 
    log = ['CSIRO GRIDDED POINTS DOWNLOAD LOG\n\n']
    for pi in points_info:
        log.append('longitude: {0:.2f}, latitude: {1:.2f}\n'.format(
            pi[2], pi[3])
        )
    print(''.join(log))
    log.append('\ndownload start: {0}\n'.format(str(datetime.now())))

    # download stations
    print('downloading CSIRO gridded points data... {0} files'.format(len(l_urls)))
    for u in l_urls:
        print(op.basename(u), ' ... ', end='')
        sys.stdout.flush()
        url_extract_gridded_points(u, points_info, p_dl)
        print('downloaded.')

    # join .nc files
    for p_ID in points_ID:
        p_st_dl = op.join(p_dl, '{0}'.format(p_ID))
        p_ncf = op.join(p_dl, '..', 'gridded_{0}.nc'.format(p_ID))
        join_ncs(p_st_dl, p_ncf)

    # end log
    log.append('download end:   {0}\n'.format(str(datetime.now())))
    with open(op.join(p_dl, '..', 'info_gridded_points.txt'), 'w') as fW:
        fW.writelines(log)

def download_gridded(p_store, point_list, grid_code='glob_24m',
                     time_limits=None):
    '''
    Locate and download CSIRO gridded separated points data, then stores it on netcdf format

    point_list - query, list of tuples (lon_point, lat_point)
    p_store    - path to store downloaded points
    grid_code = 'aus_4m', 'aus_10m', 'glob_24m', 'pac_4m', 'pac_10m'
    time_limits - optional tuple (time_start, time_end) with 'yyyymm' format
    '''

    # paths
    p_dl = op.join(p_store, 'gridded_points.download')  # temp storage of monthly files
    if not op.isdir(p_dl): os.makedirs(p_dl)

    # Generate URL list 
    l_urls = generate_urls('gridded', grid_code, time_limits=time_limits)

    # get stations lat,lon from first file
    with nc4.Dataset(l_urls[0], 'r') as ff:
        slons = ff['longitude'][:]
        slats = ff['latitude'][:]

    # unravel lonlat point combinations
    mg_lon, mg_lat = np.meshgrid(slons, slats)
    mg_lon = mg_lon.ravel()
    mg_lat = mg_lat.ravel()

    # find nearest points coordinates
    points_ix = []
    for lon_p, lat_p in point_list:
        points_ix.append((np.sqrt(np.square(mg_lat-lat_p)+np.square(mg_lon-lon_p))).argmin())

        # TODO calculate geodesic distance (degree)
        #geo_dist = []
        #for lon_ps, lat_ps in zip(slons, slats):
        #     geo_dist.append(geo_distance(lat_ps, lon_ps, p_lat, p_lon))
        #geo_dist = np.asarray(geo_dist)

    # join some station information in a list of tuples
    points_info = [
        (x, '{0:.3f}_{1:.3f}'.format(mg_lon[x], mg_lat[x]),
         mg_lon[x], mg_lat[x]) for x in points_ix
    ]

    # download selected stations
    download_gridded_points(l_urls, p_dl, points_info)

def download_gridded_area(p_store, lon1, lat1, lon2, lat2,
                          grid_code='glob_24m', time_limits=None):
    '''
    Download CSIRO gridded data and stores it on netcdf format

    stations are searched inside area:
    lon1, lat1 - bottom left
    lon2, lat2 - top right
    grid_code = 'aus_4m', 'aus_10m', 'glob_24m', 'pac_4m', 'pac_10m'
    p_store    - path to store downloaded stations
    time_limits - optional tuple (time_start, time_end) with 'yyyymm' format
    '''

    # paths
    p_dl = op.join(p_store, 'gridded_area.download')  # temp storage of monthly files
    p_ncf = op.join(p_store, 'gridded_area.nc')  # temp storage of monthly files
    if not op.isdir(p_dl): os.makedirs(p_dl)

    # Generate URL list 
    l_urls = generate_urls('gridded', grid_code, time_limits=time_limits)

    # find area to extract indexes
    with xr.open_dataset(l_urls[0]) as ff:
        idx1 = (np.abs(ff.longitude.values - lon1)).argmin()
        idy1 = (np.abs(ff.latitude.values - lat1)).argmin()
        idx2 = (np.abs(ff.longitude.values - lon2)).argmin()
        idy2 = (np.abs(ff.latitude.values - lat2)).argmin()

    # ini log 
    log = ['CSIRO GRIDDED AREA DOWNLOAD LOG\n\n']
    log.append('longitude: {0:.3f} - {1:.3f}\n'.format(
        ff.longitude.values[idx1], ff.longitude.values[idx2]))
    log.append('latitude: {0:.3f} - {1:.3f}\n'.format(
        ff.latitude.values[idy1], ff.latitude.values[idy2]))
    print(''.join(log))
    log.append('\ndownload start: {0}\n'.format(str(datetime.now())))

    # download stations
    print('downloading CSIRO gridded area data... {0} files'.format(len(l_urls)))
    for u in l_urls:
        print(op.basename(u), ' ... ', end='')
        sys.stdout.flush()
        url_extract_gridded_area(u, idx1, idy1, idx2, idy2, p_dl)
        print('downloaded.')

    # TODO: need to separate at 201306?
    # find gridded-data-format split position, split lists
    #pos_split =['201306' in x for x in l_urls].index(True)

    # join .nc files
    join_ncs(p_dl, p_ncf)

    # end log
    log.append('download end:   {0}\n'.format(str(datetime.now())))
    with open(op.join(p_dl, '..', 'info_gridded_area.txt'), 'w') as fW:
        fW.writelines(log)


# LON, LAT SPEC/GRIDDED INFO
def download_info_spec(p_ncfile):
    '''
    Download CSIRO spectral stations lat, long, name and
    store it on netcdf format.

    returns xarray.Dataset (station) latitude, longitude, station_name
    '''

    # Generate URL list 
    l_urls = generate_urls('spec')

    # get stations 
    print('downloading CSIRO spec stations...')
    with xr.open_dataset(l_urls[0]) as ff:

        # generate output dataset
        xds_out = xr.Dataset(
            {
                'longitude': (('station'), ff.longitude[0,:]),
                'latitude': (('station'), ff.latitude[0,:]),
                'station_name': (('station'), ff.station_name[:]),
            },
            coords = {
                'station' : ff.station,
            }
        )
    print('done.')

    # save to netcdf file
    xds_out.to_netcdf(p_ncfile, 'w')

    return xds_out

def download_info_gridded(p_ncfolder):
    '''
    Download CSIRO grids lat, lon and attrs
    grid_code:  'pac_4m', 'pac_10m', 'aus_4m', 'aus_10m', 'glob_24m'

    returns list of xarray.Dataset
    '''

    grid_codes = ['glob_24m', 'pac_4m', 'pac_10m', 'aus_4m', 'aus_10m']

    # download folder
    if not op.isdir(p_ncfolder):
        os.makedirs(p_ncfolder)

    # Generate URL list 
    l_xds_grids = []
    for gc in grid_codes:
        print('downloading gridded coordinates: {0} ... '.format(gc))
        l_urls = generate_urls('gridded', gc)
        with xr.open_dataset(l_urls[0]) as ff:

            # get mask
            hs1 = ff.hs.isel(time=0).values[:]
            mask = ~np.isnan(hs1)

            # generate output dataset
            xds_out = xr.Dataset(
                {
                    'mask': (('latitude', 'longitude'), mask),
                },
                coords = {
                    'longitude': ff.longitude,
                    'latitude': ff.latitude,
                },
                attrs = ff.attrs
            )
            l_xds_grids.append(xds_out)

            # save to netcdf file
            xds_out.to_netcdf(op.join(p_ncfolder,'{0}.nc'.format(gc)), 'w')

    print('done.')

    # Return grids
    return l_xds_grids

