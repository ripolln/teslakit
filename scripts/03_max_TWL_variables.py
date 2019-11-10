import sys
import  pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import datetime
from teslakit.io.matlab import ReadMatfile
from teslakit.custom_dateutils import DateConverter_Mat2Py
from teslakit.estela import Plot_DWTs_Probs
from teslakit.plotting.wts import Plot_Probs_WT_WT

rutin = '/Users/albacid/Projects/SERDP/results_files/'

#----------------------------------------------------------------
# load max TWL dataset
data = xr.open_dataset('/Users/albacid/Projects/SERDP/twl_POT.nc')

time_max = data['time'].values

#----------------------------------------------------------------
# load full dataset
# 3-hourly offshore values:
matfile = rutin + 'Historicos/KWA_historical_parameters_2016_sep.mat'

# Load data
data_full = ReadMatfile(matfile)
# print(data.keys())
del data_full['dates']

data_full = pd.DataFrame(data_full)
data_full['time'] = DateConverter_Mat2Py(data_full['time'])
data_full = data_full.set_index('time')
data_full['mmsl'] = data_full['mmsl']/1000.0 # to meters??

data_bmus = data_full['bmus']


#----------------------------------------------------------------
# Select from the bmus database, the instants of twl_max
bmus_twl_max = []

for t in time_max:

    bmus_twl_max.append(data_bmus.loc[t])
    # print(var_twl_max)

bmus_twl_max= pd.DataFrame({'time': time_max, 'bmus': bmus_twl_max})
bmus_twl_max = bmus_twl_max.set_index('time')

print(bmus_twl_max)

#----------------------------------------------------------------
# DWT probabilities associated to twl peaks
# plt.hist(bmus_twl_max['bmus'], bins=max(bmus_twl_max['bmus']))
# plt.show()

Plot_DWTs_Probs(bmus_twl_max['bmus'], bmus_twl_max.index, 42, p_export='/Users/albacid/Projects/SERDP/DWTprob_TWLpeaks.png')


#----------------------------------------------------------------
# DWT probabilities associated to twl peaks
Plot_Probs_WT_WT(
    MJO_phase, DWT_bmus, MJO_ncs, DWT_ncs,
    wt_colors=False, ttl='MJO Phases / DWT bmus (Historical)')
