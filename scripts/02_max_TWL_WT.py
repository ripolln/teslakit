import sys
import  pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import datetime
from teslakit.io.matlab import ReadMatfile
from teslakit.custom_dateutils import DateConverter_Mat2Py
from teslakit.plotting.estela import Plot_DWTs_Probs, Plot_DWTs_Mean_Anom
from teslakit.plotting.wts import Plot_Probs_WT_WT
from teslakit.custom_dateutils import xds_common_dates_daily as xcd_daily

rutin = '/Users/albacid/Projects/SERDP/results_files/'


#----------------------------------------------------------------
# load max TWL dataset
data_hist = xr.open_dataset('/Users/albacid/Projects/SERDP/twl_POT_offshore_historical.nc')

data_synth = xr.open_dataset('/Users/albacid/Projects/SERDP/twl_POT_offshore_synthetic.nc')



# #----------------------------------------------------------------
# # Plot histograms
# nbins = 8
#
# fig, ax = plt.subplots(3,2)
#
# bins_0 = ax[0,0].hist(data_hist['twl'].values, color='red', bins=nbins)
# ax[0,0].set_xlabel('twl peaks (m)')
# ax[0,0].set_ylabel('nº events')
#
# bins_1 = ax[1,0].hist(data_hist['area'].values, color='red', bins=nbins)
# ax[1,0].set_xlabel('area (m*h)')
# ax[1,0].set_ylabel('nº events')
#
# bins_2 = ax[2,0].hist(data_hist['duration'].values, color='red', bins=nbins)
# ax[2,0].set_xlabel('duration (h)')
# ax[2,0].set_ylabel('nº events')
#
# ax[0,1].hist(data_synth['twl'].values, bins=bins_0[1])
# ax[0,1].set_xlabel('twl peaks (m)')
# ax[0,1].set_ylabel('nº events')
#
# ax[1,1].hist(data_synth['area'].values, bins=bins_1[1])
# ax[1,1].set_xlabel('area (m*h)')
# ax[1,1].set_ylabel('nº events')
#
# ax[2,1].hist(data_synth['duration'].values, bins=bins_2[1])
# ax[2,1].set_xlabel('duration (h)')
# ax[2,1].set_ylabel('nº events')
#
# ax[0,0].legend('historical')
# ax[0,1].legend('synthetic')
# plt.tight_layout()
# plt.show()
# # #plt.savefig('/Users/albacid/Projects/SERDP/variables_POT_hist_synth_raw.png', dpi=600)
# # plt.savefig('/Users/albacid/Projects/SERDP/variables_POT_hist_synth.png', dpi=600)
# sys.exit()


# #----------------------------------------------------------------
# # DWT probabilities associated to twl peaks
# print(data_hist)
# print()
# print(data_synth)
# print()
#
# Plot_DWTs_Probs(data_hist['bmus'].values, data_hist['time'].values, 42, p_export='/Users/albacid/Projects/SERDP/DWTprob_TWLpeaks_hist.png')
# Plot_DWTs_Probs(data_synth['bmus'].values, data_synth['time'].values, 42, p_export='/Users/albacid/Projects/SERDP/DWTprob_TWLpeaks_synth.png')
#
# sys.exit()

#----------------------------------------------------------------
# plot area values for each DWT

data_hist['sorted_bmus'] = data_hist['bmus'] -1 # para esta funcion bmus tiene que ir de 0 a 41
data_hist['n_clusters'] = np.arange(42)

data_synth['sorted_bmus'] = data_synth['bmus'] -1 # para esta funcion bmus tiene que ir de 0 a 41
data_synth['n_clusters'] = np.arange(42)

#Plot_DWTs_Mean_Anom(data_hist, data_hist['area'], kind='mean', mask_land=None, p_export='/Users/albacid/Projects/SERDP/DWTareas_TWLpeaks_hist.png')
Plot_DWTs_Mean_Anom(data_synth, data_synth['area'], kind='mean', mask_land=None, p_export='/Users/albacid/Projects/SERDP/DWTareas_TWLpeaks_synth.png')

sys.exit()

#----------------------------------------------------------------
# DWT probabilities associated to MJO fases
MJO_hist = xr.open_dataset('/Users/albacid/Projects/TeslaKit_projects/databases/MJO/MJO_hist.nc')
MJO_phase = MJO_hist.phase.values[:]
MJO_hist['phase'] = MJO_hist.phase - 1

bmus_twl_max['bmus'] = bmus_twl_max['bmus'].values - 1
bmus_twl_max = bmus_twl_max.to_xarray()


# get common dates
dc = xcd_daily([bmus_twl_max, MJO_hist])
xds_DWT = bmus_twl_max.sel(time=slice(dc[0], dc[-1]))
xds_MJO = MJO_hist.sel(time=slice(dc[0], dc[-1]))

print(xds_MJO)
print()
print(xds_DWT)


MJO_ncs = 8
DWT_ncs = 42
Plot_Probs_WT_WT(
    xds_MJO['phase'].values[:], xds_DWT['bmus'].values[:], MJO_ncs, DWT_ncs,
    wt_colors=False, ttl='MJO Phases / DWT bmus of TWL max (Historical)')


