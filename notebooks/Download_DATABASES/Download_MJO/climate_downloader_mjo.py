#!/usr/bin/env python
# -*- coding: utf-8 -*-

# common
import time
import os
import os.path as op
from datetime import datetime, timedelta

# pip
import numpy as np
import netCDF4 as nc4
import xarray as xr
import requests


def download(p_ncfile, init_year=None):
    '''
    Download Madden Julian Oscillation data and stores it on netcdf format
    init_year: optional, data before init_year will be discarded ('yyyy-mm-dd')
    log: optional, show log

    returns xarray.Dataset
    xds_MJO:
        (time, ) mjo
        (time, ) phase
        (time, ) rmm1
        (time, ) rmm2
    '''

    # default parameter
    url_mjo = r'http://www.bom.gov.au/climate/mjo/graphics/rmm.74toRealtime.txt'
    tmp_file = 'tmp_mjo.txt'

    # download updated raw mjo table
    with open(tmp_file,'wb') as fW,  requests.get(url_mjo) as rR:
        fW.write(rR.content)

    # download data and mount time array
    ddata = np.genfromtxt(
        tmp_file,
        skip_header=2,
        usecols=(0,1,2,3,4,5,6),
        dtype = None,
        names = ('year','month','day','RMM1','RMM2','phase','amplitude'),
    )

    # mount dattime array
    dtimes = [datetime(d['year'], d['month'], d['day']) for d in ddata]

    # parse data to xarray.Dataset
    ds_mjo = xr.Dataset(
        {
            'mjo'   :(('time',), ddata['amplitude']),
            'phase' :(('time',), ddata['phase']),
            'rmm1'  :(('time',), ddata['RMM1']),
            'rmm2'  :(('time',), ddata['RMM2']),
        },
        {'time' : dtimes}
    )

    # cut dataset if asked
    if init_year:
        ds_mjo = ds_mjo.loc[dict(time=slice(init_year, None))]

    # save at netcdf file
    ds_mjo.to_netcdf(p_ncfile,'w')

    # remove temp file
    os.remove(tmp_file)

