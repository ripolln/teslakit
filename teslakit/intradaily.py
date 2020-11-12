#!/usr/bin/env python
# -*- coding: utf-8 -*-

# pip
import numpy as np
import xarray as xr

# tk
from .util.time_operations import npdt64todatetime

# TODO: refactor con waves.py/hydrographs

class Hydrograph(object):
    'Stores hydrograph data'

    def __init__(self):
        self.date_index = []
        self.dates = []

        self.indx_hydro = []
        self.numdays_hydro = []
        self.Hs_hydro = []
        self.Tp_hydro = []
        self.Dir_hydro = []

        self.TWL_max = []
        self.Hs_max = []

        self.MU = []
        self.TAU = []


def Calculate_Hydrographs(xds_BMUS, xds_WAVES):
    '''
    Calculates intradaily hydrographs

    xds_BMUS: (time) bmus
    xds_WAVES: (time) hs, tp, dir

    returns dictionary of Hydrograph objects for each WT
    and list of xarray.Datasets containing MU and TAU
    '''

    # solve intradaily bins
    bmus = xds_BMUS.bmus.values[:]
    time_KMA = xds_BMUS.time.values[:]

    d_bins = {}
    l_mutau_xdsets = []
    for i_wt in sorted(set(bmus)):

        # find WT indexes at KMA bmus
        indx = np.where((bmus == i_wt))[0]
        date_index = time_KMA[indx]

        # find hydrograms longer than 1 day
        diff_time = np.array([
            (b-a).astype('timedelta64[D]')/np.timedelta64(1,'D') \
            for a,b in zip(date_index, date_index[1:])])
        sep_hydro = np.where((diff_time > 1.0))[0]

        if len(sep_hydro)==0:
            bin_k = 'WT {0:02d}'.format(i_wt)
            d_bins[bin_k] = None
            print('{0} empty'.format(bin_k))
            continue

        hydro_indx = []
        hydro_indx.append(indx[0:sep_hydro[0]+1])
        for m in range(len(sep_hydro)-2):
            hydro_indx.append(
                indx[sep_hydro[m]+1:sep_hydro[m+1]+1]
            )
        hydro_indx.append([indx[sep_hydro[len(sep_hydro)-1]-1]])
        num_days = [len(x) for x in hydro_indx]

        # initialize some output lists
        Hs_hydro = []
        Tp_hydro = []
        Dir_hydro = []
        TWL_max = []
        Hs_max = []
        MU = []
        TAU = []

        # work with storms <= 4 days
        ndays_storms = 4
        for h, nd in zip(hydro_indx, num_days):

            # initialize loop output
            hs_W = []
            tp_W = []
            dir_W = []
            twl_max_v = 0
            twl_max_t = 0
            hs_max_v = 0
            mu_v = 0

            if nd <= ndays_storms:

                # start and end dates for hydrograph
                p1 = npdt64todatetime(time_KMA[h[0]])
                p2 = npdt64todatetime(time_KMA[h[-1]]).replace(hour=23)

                # get waves conditions for hydrograph
                xds_W = xds_WAVES.sel(time = slice(p1, p2))
                hs_W = xds_W.Hs.values[:]
                tp_W = xds_W.Tp.values[:]
                dir_W = xds_W.Dir.values[:]

                #TODO: añadido por AlbaC
                if 'runup' in xds_WAVES:
                    # calculate runup max and normalize
                    twl_temp = xds_W.runup.values[:]

                else:
                    # calculate TWL max and normalize
                    twl_temp = 0.1*(hs_W**0.5)*tp_W


                dt = 1.0/len(twl_temp)
                t_norm = np.arange(dt,1+dt,dt)

                # calculate TAU
                i_twlmax = np.where(twl_temp == np.amax(twl_temp))[0]
                twl_max_v =  twl_temp[i_twlmax][0]
                twl_max_t =  t_norm[i_twlmax][0]

                # calculate MU
                mu_v = np.trapz(np.divide(twl_temp, twl_max_v), t_norm)

                # calculate max Hs and Tp
                hs_max_v = hs_W[np.where(hs_W == np.amax(hs_W))[0]][0]
                tp_max_v = tp_W[np.where(tp_W == np.amax(tp_W))[0]][0]

            # store values
            Hs_hydro.append(hs_W)
            Tp_hydro.append(tp_W)
            Dir_hydro.append(dir_W)
            TWL_max.append(twl_max_v)
            Hs_max.append(hs_max_v)
            MU.append(mu_v)
            TAU.append(twl_max_t)

        # bin key and hydrograph storage
        bin_k = 'bin{0:02d}'.format(i_wt)
        bin_hy = Hydrograph()

        bin_hy.date_index = indx
        bin_hy.dates = date_index
        bin_hy.indx_hydro = hydro_indx
        bin_hy.numdays_hydro = num_days
        bin_hy.Hs_hydro = Hs_hydro
        bin_hy.Tp_hydro = Tp_hydro
        bin_hy.Dir_hydro = Dir_hydro
        bin_hy.TWL_max = TWL_max
        bin_hy.Hs_max = Hs_max
        bin_hy.MU = MU
        bin_hy.TAU = TAU

        # store at dictionary
        d_bins[bin_k] = bin_hy
        #print('{0} calculated'.format(bin_k))

        # store mu, tau xarray.Dataset
        xds_hg_wt = xr.Dataset(
            {
                'MU':(('time',), np.array(MU)),
                'TAU':(('time',), np.array(TAU)),
                'hs_max':(('time',), np.array(Hs_max)),
                'twl_max':(('time',), np.array(TWL_max)),
            },
            attrs={'WT':i_wt+1}
        )
        l_mutau_xdsets.append(xds_hg_wt)

    return d_bins, l_mutau_xdsets


def Calculate_Hydrographs_runup_hist(xds_BMUS, xds_WAVES, xds_RUNUP):
    '''
    Calculates intradaily hydrographs

    xds_BMUS: (time) bmus --> daily data
    xds_WAVES: (time) hs, tp, dir --> hourly data

    returns dictionary of Hydrograph objects for each WT
    and list of xarray.Datasets containing MU and TAU
    '''

    # solve intradaily bins
    bmus = xds_BMUS.bmus.values[:]
    time_KMA = xds_BMUS.time.values[:]

    dif_bmus = np.diff(bmus)
    dif_bmus = np.nonzero(dif_bmus)[0]
    dif_bmus = np.append(dif_bmus, len(bmus))

    # initialize some output lists
    STORM_wt = []
    TWL_max = []
    Hs_max = []
    Tp_max = []
    Dir_max = []
    MU = []
    MU_2 = []
    TAU = []
    T_ini = []
    T_runup_max = []


    # separate in storms
    ind_ini = 0
    for d in dif_bmus:
        ind_fin = d+1

        storm_WT = bmus[ind_ini:ind_fin]
        time_WT = time_KMA[ind_ini:ind_fin]

        #--------------------------------------------------------------------
        # # find hydrograms longer than 1 day
        # if len(storm_WT) == 1:
        #     twl_max_v = np.nan
        #     twl_max_t = np.nan
        #     hs_max_v = np.nan
        #     tp_max_v = np.nan
        #     dir_max_v = np.nan
        #     mu_v = np.nan
        #     mu_v_2 = np.nan
        #     time_ini = time_WT[0]
        #     time_fin = time_WT[0]
        #     time_runup_max = time_WT[0]


        #else:

        # start and end dates for hydrograph
        p1 = npdt64todatetime(time_WT[0])
        p2 = npdt64todatetime(time_WT[-1]).replace(hour=23)

        # get waves conditions for hydrograph
        xds_W = xds_WAVES.sel(time = slice(p1, p2))


        # calculate TWL max and normalize
        xds_R = xds_RUNUP.sel(time = slice(p1, p2))
        twl_temp = xds_R.values[:]


        dt = 1.0/len(twl_temp)
        t_norm = np.arange(dt,1+dt,dt)


        # calculate TAU
        i_twlmax = np.where(twl_temp == np.amax(twl_temp))[0]
        twl_max_v =  twl_temp[i_twlmax][0]
        twl_max_t =  t_norm[i_twlmax][0]

        time_runup_max = xds_R.time.isel(time=i_twlmax[0]).values


        # calculate MU
        mu_v = np.trapz(np.divide(twl_temp, twl_max_v), t_norm) # area debajo de la curva

        mu_v_2 = np.trapz(np.divide([twl_temp[0], twl_temp[i_twlmax], twl_temp[-1]], twl_max_v),
                        [t_norm[0], t_norm[i_twlmax], t_norm[-1]])[0] # area debajo de los 3 puntos


        # calculate associated Hs, Tp and Dir
        hs_max_v = xds_W.Hs.isel(time=i_twlmax[0]).values
        tp_max_v = xds_W.Tp.isel(time=i_twlmax[0]).values
        dir_max_v = xds_W.Dir.isel(time=i_twlmax[0]).values


        #--------------------------------------------------------------------

        # store values
        TWL_max.append(twl_max_v)
        Hs_max.append(hs_max_v)
        Tp_max.append(tp_max_v)
        Dir_max.append(dir_max_v)
        MU.append(mu_v)
        MU_2.append(mu_v_2)
        TAU.append(twl_max_t)
        T_runup_max.append(time_runup_max)
        T_ini.append(p1)
        STORM_wt.append(storm_WT[0])

        # update index
        ind_ini = ind_fin


    # store mu, tau xarray.Dataset
    xds_hg_wt = xr.Dataset(
        {
            'MU':(('time',), np.array(MU)),
            'MU2':(('time',), np.array(MU_2)),
            'TAU':(('time',), np.array(TAU)),
            'Hs':(('time',), np.array(Hs_max)),
            'Tp':(('time',), np.array(Tp_max)),
            'Dir':(('time',), np.array(Dir_max)),
            'Runup_max':(('time',), np.array(TWL_max)),
            'time_runup_max':(('time',),  np.array(T_runup_max).astype('datetime64[ns]')),
            'WT':(('time',), np.array(STORM_wt)),
        },
        coords={
                'time': np.array(T_ini).astype('datetime64[ns]'),
                },
    )


    return xds_hg_wt
