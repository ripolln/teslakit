
import os.path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import sys
from teslakit.io.matlab import ReadMatfile
from datetime import datetime, timedelta


#--------------------------------------------------------------------------------
def datematlab2datetime(datenum_matlab):
    'Return python datetime for matlab datenum. Transform and adjust from matlab.'

    d_new = list(np.ones(len(datenum_matlab),) * np.nan)

    for i, d in enumerate(datenum_matlab):
        d_temp = datetime.fromordinal(int(d)) + \
    timedelta(days=float(d % 1)) - \
        timedelta(days=366) + timedelta(microseconds=0)

        d_new[i] = d_temp

    return d_new


# #------------------------------------------------------------------------
# # 1) read TG in txt, save netcdf
# rutin = '/Users/albacid/Projects/TeslaKit_projects/databases/Tide/Guam/'
# file = '1630000.19960101.20161231.x.xx.HHt.OPR.txt'
#
# data = pd.read_csv(
#     os.path.join(rutin, file),
#     names = ['i', 'date', 'observed', 'AT', 'NTR', 'sigma'],
#     sep='[s+\t]', engine='python'
# )
#
# data = data.drop(columns='i')
#
# time = pd.to_datetime(data['date'])
# data['date'] = time
#
# data = data.set_index('date')
#
# ds = data.to_xarray()
#
# ds.to_netcdf(os.path.join(rutin, 'GUAM_TG_1996_2016.nc'))
#
# print(ds)
# sys.exit()


# #------------------------------------------------------------------------
# # 2) Read TG in mat file, save to netcdf
# rutin = '/Users/albacid/Projects/TeslaKit_projects/sites/z_GUAM_Laura/TIDE'
# file = 'Mareografo_Guam.mat'
# TG = ReadMatfile(os.path.join(rutin, file))
#
# time = datematlab2datetime(TG['DATES'])
#
# xds_out = xr.Dataset({ }, coords = {'time': time},)
# for k in TG.keys():
#
#     xds_out[k] =(('time',), TG[k])
#
# xds_out = xds_out.drop('DATES')
# print(xds_out.isel(time=slice(0,2)))
#
# xds_out = xds_out.resample(time='H').mean()
# print(xds_out.NTR.isel(time=slice(0,2)))
# sys.exit()
#
# xds_out.to_netcdf(os.path.join(rutin,'GUAM_TG_1989_2018.nc'))
# sys.exit()


#------------------------------------------------------------------------
# 3) plot data
rutin = '/Users/albacid/Projects/TeslaKit_projects'
TG_sigma = xr.open_dataset(os.path.join(rutin, 'databases', 'Tide', 'Guam','GUAM_TG_1996_2016.nc'))
#TG = xr.open_dataset(os.path.join(rutin, 'databases', 'Tide', 'Guam','GUAM_TG_1989_2018.nc'))

print(TG_sigma)


# plot figure
fig, axs = plt.subplots(figsize=(12,9))
plt.plot(TG_sigma.date, TG_sigma.observed, '-r', linewidth = 0.04)
plt.plot(TG_sigma.date, TG_sigma.AT, '-b', linewidth = 0.04)
plt.plot(TG_sigma.date, TG_sigma.NTR, '-k', linewidth = 0.04)
plt.plot(TG_sigma.date, TG_sigma.sigma, '-m', linewidth = 0.04)

#plt.xlim(TG_sigma.date[0], TG_sigma.date[-1])
#plt.title('Astronomical tide')
plt.xlabel('time')
#plt.ylabel('tide (m)')
plt.show()

