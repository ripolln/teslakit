import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import sys


# # ------------------------------------------
# nc = '/Users/albacid/Projects/TeslaKit_projects/sites/GUAM/TCs/Allstorms.ibtracs_wmo.v03r10.nc'
nc = '/Users/albacid/Projects/TeslaKit_projects/sites/GUAM/TCs/TCs_hist_r1.nc'
# nc = '/Users/albacid/Projects/TeslaKit_projects/sites/KWAJALEIN/TCs/TCs_hist_r1.nc'


TCs = xr.open_dataset(nc)

lon = TCs.pres_wmo.values[:]

lon_end = np.empty(len(TCs.storm))*np.nan

for s in TCs.storm:

    lon_s = lon[s]

    ind_nonan = np.where(~np.isnan((lon_s)))
    lon_s = lon_s[ind_nonan]

    lon_end[s] = np.min(lon_s)

ind = np.where(lon_end==np.min(lon_end))
print(lon_end[ind])
print(ind[0])
print(TCs.pres_wmo.isel(storm=ind[0][0]))
print('------')

nc = '/Users/albacid/Projects/TeslaKit_projects/sites/GUAM/TCs/TCs_hist_r1_params.nc'
TCs = xr.open_dataset(nc)
print(TCs.pressure_min.isel(storm=ind[0][0]))

print(TCs.category.isel(storm=ind[0][0]))
sys.exit()


# # ------------------------------------------
# # nc = '/Users/albacid/Projects/TeslaKit_projects/sites/KWAJALEIN/TCs/TCs_hist_r1.nc'
# nc = '/Users/albacid/Projects/TeslaKit_projects/sites/GUAM/TCs/Allstorms.ibtracs_wmo.v03r10.nc'
#
# TCs = xr.open_dataset(nc)
# lon = TCs.lon_wmo.values[:]
#
# lon_end = np.empty(len(TCs.storm))*np.nan
#
# for s in TCs.storm:
#
#     lon_s = lon[s]
#
#     ind_nonan = np.where(~np.isnan((lon_s)))
#     lon_s = lon_s[ind_nonan]
#
#     lon_end[s] = lon_s[-1]
#
#
# plt.plot(lon_end, '.')
# plt.show()
# sys.exit()

# # ------------------------------------------
# nc = '/Users/albacid/Projects/TeslaKit_projects/sites/GUAM/TCs/Allstorms.ibtracs_wmo.v03r10.nc'
# TCs = xr.open_dataset(nc)
# lon = TCs.lon_wmo.values[:]
#
# storm = 2 # datos diarios
# print(TCs['time_wmo'][storm][:10].values)
# print()
#
# storm = 10 # datos 6-h
# print(TCs['time_wmo'][storm][:10].values)
# print()
#
# storm = 1000 # datos irregulares•
# print(TCs['time_wmo'][storm][:10].values)

