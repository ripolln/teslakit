#!/usr/bin/env python
# -*- coding: utf-8 -*-

# common
import sys
import os
import os.path as op
import urllib.request
import requests
import gzip
import shutil
import time
from bs4 import BeautifulSoup

# pip
import numpy as np
import xarray as xr


def download_ibtracs_all(p_ncfile):
    '''
    Download storms (Allstorms.ibtracs_all) netcdf from NOAA

    ftp file name: Allstorms.ibtracs_all.v03rXX.nc.
    ftp wind velocity in knots: x1.82 km/h
    '''

    # default parameters
    # TODO: update version
    ftp_down = 'ftp://eclipse.ncdc.noaa.gov/pub/ibtracs/v03r10/wmo/netcdf/'
    fil_down = 'Allstorms.ibtracs_wmo.v03r10.nc.gz'

    p_down = op.dirname(p_ncfile)
    p_gz = op.join(p_down, '{0}'.format(fil_down))
    p_temp = op.join(p_down, 'tempdl.nc')  # temporal file

    # download gz
    ftp_wmo = '{0}{1}'.format(ftp_down, fil_down)
    if not op.isfile(p_gz):
        urllib.request.urlretrieve(ftp_wmo, p_gz)

    # decompress .gz file
    with gzip.open(p_gz, 'rb') as f_in, open(p_temp, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

    # load and return xarray.Dataset
    xds_wmo = xr.open_dataset(p_temp).copy()

    # set lon to 0-360
    lon_wmo = xds_wmo.lon_wmo.values[:]
    lon_wmo[np.where(lon_wmo<0)] = lon_wmo[np.where(lon_wmo<0)]+360
    xds_wmo['lon_wmo'].values[:] = lon_wmo

    # store changes
    xds_wmo.to_netcdf(p_ncfile, 'w')

    # remove temp files
    os.remove(p_temp)
    os.remove(p_gz)

def download_sst_v5(p_store, overwrite=False):
    '''
    Download NOAA SST data and stores it on netcdf format

    p_store    - path to store downloaded data
    overwrite  - True for overwriting previously downloaded files
    '''
    # prepare download folder
    p_dl = op.join(p_store, 'sst_v5.download')  # temp storage of monthly files
    if not op.isdir(p_dl): os.makedirs(p_dl)

    # default parameters
    url = 'https://www1.ncdc.noaa.gov/pub/data/cmb/ersst/v5/netcdf/'
    code = 'ersst.v5'

    # get files at public folder
    page = requests.get(url)
    soup = BeautifulSoup(page.text, 'html.parser')
    fs = ['{0}{1}'.format(url, n.get('href')) for n in soup.find_all('a') if n.get('href').startswith(code)]

    # download files
    print('downloading files ')
    for remote in fs:

        fn = op.basename(remote)
        local = op.join(p_dl, fn)

        # skip if overwrite off
        if overwrite == False and op.isfile(local):
            continue

        with requests.get(remote) as fr, open(local, 'wb')as fl:
            fl.write(fr.content)
        print(fn)

    print()

    # join files
    fjs = [op.join(p_dl, x) for x in sorted(os.listdir(p_dl)) if x.endswith('.nc')]
    ts = [d.split('.')[-2] for d in fjs]
    ts = [np.datetime64('{0}-{1}-{2}'.format(x[:4], x[4:], '01')) for x in ts]

    # get coordinates and attributes from downloaded files
    tmp_0 = xr.open_dataset(fjs[0]).squeeze(drop=True)
    tmp_1 = xr.open_dataset(fjs[-1]).squeeze(drop=True)

    # initialize output file
    out = xr.Dataset(coords=tmp_0.coords)
    out['time'] = ts  # update time attribute
    out.attrs = tmp_1.attrs
    out.attrs['id'] = code
    out.attrs['time_coverage_start'] = tmp_0.attrs['time_coverage_start']

    # initialize dataset variables
    sh_vs = tmp_0['sst'].shape + (len(ts),)
    dvs = {
        'sst': np.zeros(sh_vs)*np.nan,
        'ssta': np.zeros(sh_vs)*np.nan,
    }

    # read variables
    for c, k in enumerate(fjs):
        xk = xr.open_dataset(k)

        for vn in ['sst', 'ssta']:
            dvs[vn][:,:,c] = xk[vn].values[:]

    # assign variables
    for vn in ['sst', 'ssta']:
        out[vn] = (('lat', 'lon', 'time',), dvs[vn])

    # store output
    mf = op.join(p_store, '{0}.nc'.format(code))
    out.to_netcdf(mf)
    print('Complete. Data stored at: {0}'.format(mf))

    return out

