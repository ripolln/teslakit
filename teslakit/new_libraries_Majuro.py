#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import xarray as xr
from datetime import datetime, timedelta

# tk
from .util.time_operations import get_years_months_days, npdt64todatetime, \
fast_reindex_hourly

# hide numpy warnings
np.warnings.filterwarnings('ignore')

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import xarray as xr
from datetime import datetime, timedelta

# tk
from .util.time_operations import get_years_months_days, npdt64todatetime, \
fast_reindex_hourly

# hide numpy warnings
np.warnings.filterwarnings('ignore')

def GetDistribution_ws(xds_wps, swell_sectors):
    '''
    Separates wave partitions (0-5) into families.
    Default: sea, swl1, swl2

    xds_wps (waves partitionss):
        xarray.Dataset (time,), phs, pspr, pwfrac... {0-n partitions}

    sectors: list of degrees to cut wave energy [(a1, a2), (a2, a3), (a3, a1)]

    returns
        xarray.Dataset (time,), fam_V, {fam: sea,swell_1,swell2. V: Hs,Tp,Dir}
    '''

# concatenate energy groups
    sea_Hs= xds_wps.isel(part=0).hs.values
    sea_Tp=xds_wps.isel(part=0).tp.values
    sea_Dir = xds_wps.isel(part=0).dpm.values
    time= xds_wps.time.values
    
    cat_hs=np.full([len(xds_wps.time),len(xds_wps.part)-1],np.nan)
    for a in range(len(xds_wps.part)-1):
        cat_hs[:,a]=xds_wps.isel(part=a+1).hs.values
    
    cat_tp=np.full([len(xds_wps.time),len(xds_wps.part)-1],np.nan)
    for a in range(len(xds_wps.part)-1):
        cat_tp[:,a]=xds_wps.isel(part=a+1).tp.values
    
    cat_dir=np.full([len(xds_wps.time),len(xds_wps.part)-1],np.nan)
    for a in range(len(xds_wps.part)-1):
        cat_dir[:,a]=xds_wps.isel(part=a+1).dpm.values

    # prepare output array
    xds_parts = xr.Dataset({
        'sea_Hs':('time',sea_Hs),
        'sea_Tp':('time',sea_Tp),
        'sea_Dir':('time',sea_Dir)
    },
        coords = {'time':time}
    )

    # solve sectors
    c = 1
    for s_ini, s_end in swell_sectors:
        if s_ini < s_end:
            p_sw = np.where((cat_dir <= s_end) & (cat_dir > s_ini))
        else:
            p_sw = np.where((cat_dir <= s_end) | (cat_dir > s_ini))

        # get data inside sector
        sect_dir = np.zeros(cat_dir.shape)*np.nan
        sect_hs = np.zeros(cat_dir.shape)*np.nan
        sect_tp = np.zeros(cat_dir.shape)*np.nan

        sect_dir[p_sw] = cat_dir[p_sw]
        sect_hs[p_sw] = cat_hs[p_sw]
        sect_tp[p_sw] = cat_tp[p_sw]

        # calculate swell Hs, Tp, Dir
        swell_Hs = np.sqrt(np.nansum(np.power(sect_hs,2), axis=1))

        swell_Tp = np.sqrt(
            np.nansum(np.power(sect_hs,2), axis=1) /
            np.nansum(np.power(sect_hs,2)/np.power(sect_tp,2), axis=1)
        )
        swell_Dir = np.arctan2(
            np.nansum(np.power(sect_hs,2) * sect_tp * np.sin(sect_dir*np.pi/180), axis=1),
            np.nansum(np.power(sect_hs,2) * sect_tp * np.cos(sect_dir*np.pi/180), axis=1)
        )

        # dir correction and denormalization 
        swell_Dir[np.where((swell_Dir<0))] = swell_Dir[np.where((swell_Dir<0))]+2*np.pi
        swell_Dir = swell_Dir*180/np.pi

        # dont do arctan2 if there is only one dir
        i_onedir = np.where(
            (np.count_nonzero(~np.isnan(sect_dir),axis=1)==1)
        )[0]
        swell_Dir[i_onedir] = np.nanmin(sect_dir[i_onedir], axis=1)

        # out of bound dir correction
        swell_Dir[np.where((swell_Dir>360))] = swell_Dir[np.where((swell_Dir>360))]-360
        swell_Dir[np.where((swell_Dir<0))] = swell_Dir[np.where((swell_Dir<0))]+360


        # fix swell all-nans to 0s nansum behaviour
        p_fix = np.where(swell_Hs==0)
        swell_Hs[p_fix] = np.nan
        swell_Dir[p_fix] = np.nan

        # append data to partitons dataset
        xds_parts['swell_{0}_Hs'.format(c)] = ('time', swell_Hs)
        xds_parts['swell_{0}_Tp'.format(c)] = ('time', swell_Tp)
        xds_parts['swell_{0}_Dir'.format(c)] = ('time', swell_Dir)
        c+=1

    return xds_parts