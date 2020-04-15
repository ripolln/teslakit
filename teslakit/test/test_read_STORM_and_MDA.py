

#!/usr/bin/env python
# -*- coding: utf-8 -*-

# common
import os
import os.path as op
import glob
import pandas as pd

# pip
import xarray as xr
import numpy as np

# DEV: override installed teslakit
import sys
sys.path.insert(0, op.join(os.path.abspath(''), '..', '..', '..'))

# teslakit
from teslakit.database import Database
from teslakit.storms import Extract_Circle_STORM
from teslakit.mda import MaxDiss_Simplified_NoThreshold
from teslakit.plotting.storms import Plot_TCs_STORMSbase_Tracks, Plot_TCs_Params_STORM_MDAvsSIM, \
    Plot_TCs_Params_HISTvsSIM, Plot_TCs_Params_HISTvsSIM_histogram

# %%
# --------------------------------------
# Teslakit database

p_data = r'/Users/anacrueda/Documents/Proyectos/TESLA/data' # ROI/MAJURO
p_data_temp = r'/Users/anacrueda/Documents/Data/STORMs/data'

db = Database(p_data)

# set site
db.SetSite('MAJURO')

# %%
# --------------------------------------
# load data and set parameters

# choose the basin

basin = 'WP'

path = os.path.join(p_data_temp, "STORM_DATA_IBTRACS_{}*.txt".format(basin))
file_list = sorted(glob.glob(path))

tc_list = []
cont = 0
for file in file_list[0:1]:
    df = pd.read_table(file,
                       sep = ',',
                       header = None,
                       names=('year', 'month', 'storm', 'TimeStep', 'Basin', 'lat', 'lon', 'pressure_min','windspeed_max','rmax','cat','land','dist_land'),
                       dtype={'year': np.float64,
                              'month': np.float64,
                              'storm': np.float64,
                              'TimeStep': np.float64,
                              'basin': np.int64,
                              'lat': np.float64,
                              'lon': np.float64,
                              'pressure_min': np.float64,
                              'windspeed_max': np.float64,
                              'rmax': np.float64,
                            'cat': np.int64,
                            'land': np.float64,
                            'dist_land': np.float64})

    #modify storm column to know how many ciclones per simulation
    for i in df['year'].unique():
        if i == 0: continue
        late_storm = np.nanmax(df.loc[df.year == i-1, ['storm']])
            #print(late_storm)
        df.loc[df.year == i, ['storm']] = df.loc[df.year == i, ['storm']] + late_storm

        if cont > 0:
            last_storm = storms[-1]
        #print('last_storm')
        #print(last_storm)
            df['storm'] = df['storm']+last_storm
    df.set_index('storm', inplace=True)
    # store DataFrame in list
    tc_list.append(df)
    storms = df.index.tolist()
    cont = cont + 1

    # see pd.concat documentation for more info
    tc_list = pd.concat(tc_list)
tc_storms = tc_list.to_xarray()

print(tc_storms)

##
# wave point longitude and latitude ROI
#pnt_lon = 167.5
#pnt_lat = 9.75

# wave point longitude and latitude MAJURO
pnt_lon = 171.25
pnt_lat = 7.10

# radius for TCs selection (º)
r1 = 14
r2 = 4

# Get STORMS TCs at the study area
# dictionary with needed variable names
d_vns = {
    'longitude':'lon',
    'latitude':'lat',
    'time': 'TimeStep',
    'pressure':'pressure_min',
    'radius':'rmax',
    'mwinds': 'windspeed_max',
}

# Extract STORMS TCs inside r2
TCs_r2_sim_tracks, TCs_r2_sim_params = Extract_Circle_STORM(tc_storms, pnt_lon, pnt_lat, r2, d_vns)
print(TCs_r2_sim_params)

# Store STORMS TCs parameters
db.Save_TCs_r2_sim_params(TCs_r2_sim_params)


# Plot storm tracks world map (requires basemap)

lon1, lon2 = 90, 270
lat1, lat2 = -20, 70

Plot_TCs_STORMSbase_Tracks(
    TCs_r2_sim_tracks,
    lon1, lon2, lat1, lat2,
    pnt_lon, pnt_lat,  r2,
)


## check STORMs TCs vs. historical

db.Load_TCs_r2_sim_params()

# Historical vs STORM parameters:
_, TCs_r2_hist_params = db.Load_TCs_r2_hist()  # historical TCs parameters inside radius 2

# scatter plot
Plot_TCs_Params_HISTvsSIM(TCs_r2_hist_params, TCs_r2_sim_params)

# histogram
Plot_TCs_Params_HISTvsSIM_histogram(TCs_r2_hist_params, TCs_r2_sim_params)


## MDA selection

# --------------------------------------
# MaxDiss classification

# MDA number of cases
num_sel_mda = 400

# get simulated parameters
pmean_s = TCs_r2_sim_params.pressure_mean.values[:]
pmin_s = TCs_r2_sim_params.pressure_min.values[:]
gamma_s = TCs_r2_sim_params.gamma.values[:]
delta_s = TCs_r2_sim_params.delta.values[:]
vmean_s = TCs_r2_sim_params.velocity_mean.values[:]
rmax_s = TCs_r2_sim_params.mean_radius.values[:]
winds_s = TCs_r2_sim_params.winds_mean.values[:]

# subset, scalar and directional indexes
data_mda = np.column_stack((pmean_s, pmin_s, vmean_s, rmax_s, winds_s,delta_s, gamma_s))
ix_scalar = [0,1,2,3,4]
ix_directional = [5,6]

# MDA
centroids = MaxDiss_Simplified_NoThreshold(
    data_mda, num_sel_mda, ix_scalar, ix_directional
)


# store MDA storms - parameters
TCs_r2_MDA_params = xr.Dataset(
    {
        'pressure_mean':(('storm'), centroids[:,0]),
        'pressure_min':(('storm'), centroids[:,1]),
        'velocity_mean':(('storm'), centroids[:,2]),
        'mean_radius':(('storm'), centroids[:,3]),
        'winds_mean':(('storm'), centroids[:,4]),
        'delta':(('storm'), centroids[:,5]),
        'gamma':(('storm'), centroids[:,6]),
    },
    coords = {
        'storm':(('storm'), np.arange(num_sel_mda))
    },
)

print(TCs_r2_MDA_params)

db.Save_TCs_r2_mda_params(TCs_r2_MDA_params)

#  Simulated vs MDA selection: scatter plot parameters
Plot_TCs_Params_STORM_MDAvsSIM(TCs_r2_MDA_params, TCs_r2_sim_params)

