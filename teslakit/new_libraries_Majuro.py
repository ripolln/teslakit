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


def Simulate_TCs_NoEnter(self, xds_DWT, WVS_sims, xds_TCs_params,
                     xds_TCs_simulation, prob_change_TCs, MU_WT, TAU_WT):
        '''
        Climate Emulator DWTs TCs simulation

        xds_DWT             - xarray.Dataset, vars: evbmus_sims (time,)
        WVS_sim             - xarray.Dataset, output from Simulate_Waves()

        xds_TCs_params      - xr.Dataset. vars(storm): pressure_min
        xds_TCs_simulation  - xr.Dataset. vars(storm): mu, hs, ss, tp, dir
        prob_change_TCs     - cumulative probabilities of TC category change
        MU_WT, TAU_WT       - intradaily hidrographs for each WT
        '''

        # max. storm waves and KMA
        xds_KMA_MS = self.KMA_MS

        # vars needed
        dwt_bmus_sim = xds_DWT.evbmus_sims.values[:]
        dwt_time_sim = xds_DWT.time.values[:]
        n_clusters = len(xds_KMA_MS.n_clusters)

        # iterate waves simulations 
        ls_tcs_sim = []
        ls_wvs_upd = []
        for i_sim in WVS_sims.n_sim:
            wvs_s = WVS_sims.sel(n_sim=i_sim)

            # generate TCs
            tcs_sim, wvs_upd_sim = self.GenerateTCs_NoEnter(
                n_clusters, dwt_bmus_sim, dwt_time_sim,
                xds_TCs_params, xds_TCs_simulation, prob_change_TCs, MU_WT, TAU_WT,
                wvs_s
            )
            ls_tcs_sim.append(tcs_sim)
            ls_wvs_upd.append(wvs_upd_sim)

        # concatenate simulations 
        TCs_sim = xr.concat(ls_tcs_sim, 'n_sim')
        WVS_upd = xr.concat(ls_wvs_upd, 'n_sim')

        return TCs_sim, WVS_upd	

def GenerateTCs_NoEnter(self, n_clusters, DWT, DWT_time,
                    TCs_params, TCs_simulation, prob_TCs, MU_WT, TAU_WT,
                    xds_wvs_sim):
        '''
        Climate Emulator DWTs TCs simulation

        n_clusters      - KMA number of clusters
        DWT             - np.array with DWT bmus sim series (dims: time,)

        TCs_params      - xr.Dataset. vars(storm): pressure_min
        TCs_simulation  - xr.Dataset. vars(storm): mu, hs, ss, tp, dir
        prob_TCs        - cumulative probabilities of TC category change
        MU_WT, TAU_WT   - intradaily hidrographs for each WT
        xds_wvs_sim     - xr.Dataset, waves simulated without TCs (for updating)

        returns xarray.Datasets with updated Waves and simulated TCs data
            vars waves:
                *fam*_*vn* (fam: sea, swell_1, swell_2 ..., vn: Hs, Tp, Dir),
            vars TCS:
                mu, tau, ss
            dims: storm
        '''

        # wave family to modify
        mod_fam = 'sea'  # TODO input parameter

        # waves families - variables (sorted for simulation output)
        wvs_fams = self.fams
        wvs_fams_vars = [
            ('{0}_{1}'.format(f,vn)) for f in wvs_fams for vn in['Hs', 'Tp','Dir']
            ]

        # simulate one value for each storm 
        dwt_df = np.diff(DWT)
        dwt_df[-1] = 1  # ensure last day storm
        ix_ch = np.where((dwt_df != 0))[0]+1
        ix_ch = np.insert(ix_ch, 0, 0)  # get first day storm
        DWT_sim = DWT[ix_ch]
        DWT_time_sim = DWT_time[ix_ch]

        # get simulated waves for updating
        sim_wvs = np.column_stack([
            xds_wvs_sim[vn].values[:] for vn in wvs_fams_vars
        ])

        # new progress bar 
        pbar = tqdm(
            total=len(DWT_sim),
            desc = 'C.E: Sim. TCs  '
        )

        # Simulate TCs (mu, ss, tau)
        sims_out = np.zeros((len(DWT_sim), 3))
        c = 0
        while c < len(DWT_sim):
            WT = int(DWT_sim[c])
            iwt = WT - 1

            # KMA Weather Types tcs generation
            if WT <= n_clusters:

                # get random MU,TAU from current WT
                # TODO: random excluyente?
                ri = randint(len(MU_WT[iwt]))
                mu_s = MU_WT[iwt][ri]
                tau_s = TAU_WT[iwt][ri]
                ss_s = 0

            # TCs Weather Types waves generation
            else:
                
                # TC does not enter. random mu_s, 0.5 tau_s, 0 ss_s
                all_MUs = np.concatenate(MU_WT)
                ri = randint(len(all_MUs))
                # TODO: check mu 0s, set nans (?)

                mu_s = all_MUs[ri]
                tau_s = 0.5
                ss_s = 0

            sim_row = np.array([mu_s, tau_s, ss_s])

            # no nans or values < 0 stored 
            if ~np.isnan(sim_row).any() and len(np.where(sim_row<0)[0])==0:
                sims_out[c] = sim_row
                c+=1

                # progress bar
                pbar.update(1)

        pbar.close()

        # update waves simulation
        xds_WVS_sim_updated = xr.Dataset(
            {
                'DWT': (('time',), DWT_sim),
            },
            coords = {'time': DWT_time_sim}
        )
        for c, vn in enumerate(wvs_fams_vars):
            xds_WVS_sim_updated[vn] = (('time',), sim_wvs[:,c])

        # generated TCs 
        xds_TCs_sim = xr.Dataset(
            {
                'mu':  (('time',), sims_out[:,0]),
                'tau': (('time',), sims_out[:,1]),
                'ss':  (('time',), sims_out[:,2]),

                'DWT': (('time',), DWT_sim),
            },
            coords = {'time': DWT_time_sim}
        )

        return xds_TCs_sim, xds_WVS_sim_updated