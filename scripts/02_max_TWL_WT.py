import sys
import  pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from _datetime import datetime, timedelta
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
#
# #plt.show()
#
# #plt.savefig('/Users/albacid/Projects/SERDP/variables_POT_hist_synth_raw.png', dpi=600)
# plt.savefig('/Users/albacid/Projects/SERDP/variables_POT_hist_synth.png', dpi=600)
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


# #--------------------------------
# # plot area values for each DWT
# from teslakit.util.operations import GetBestRowsCols
# import matplotlib.gridspec as gridspec
# from teslakit.plotting.config import _faspect, _fsize, _fdpi
# from teslakit.kma import ClusterProbabilities
# from teslakit.plotting.wts import axplot_WT_Probs, axplot_WT_Hist
# from teslakit.plotting.custom_colors import colors_dwt
#
# def Plot_DWTs_Mean_Anom_modif(xds_KMA, xds_var, kind='mean', p_export=None):
#     '''
#     Plot Daily Weather Types (bmus mean)
#     kind - mean/anom
#     '''
#
#     bmus = xds_KMA['sorted_bmus'].values[:]
#     n_clusters = len(xds_KMA.n_clusters.values[:])
#
#     # plot figure
#     fig, ax = plt.subplots()
#
#     values = []
#     for ic in range(n_clusters):
#
#         if kind=='mean':
#             # data mean
#             it = np.where(bmus==ic)[0][:]
#             c_plot = xds_var.isel(time=it).mean(dim='time')
#
#         elif kind=='anom':
#             # data anomally
#             it = np.where(bmus==ic)[0][:]
#             t_mean = xds_var.mean(dim='time')
#             c_mean = xds_var.isel(time=it).mean(dim='time')
#             c_plot = c_mean - t_mean
#
#         values = np.append(values, c_plot.values)
#
#     values = np.reshape(values,(7,6))
#     values = np.flipud(values)
#
#     plt.pcolor(values)
#     plt.colorbar()
#     plt.clim(0,8)
#     ax.tick_params(axis='both', which='both', length=0)
#     ax.set_xticklabels('')
#     ax.set_yticklabels('')
#
#     plt.title('mean Area below TWL max (m*h)')
#
#     # show / export
#     if not p_export:
#         plt.show()
#     else:
#         fig.savefig(p_export, dpi=_fdpi)
#         plt.close()
#
# #--------------------------------
#
#
# data_hist['sorted_bmus'] = data_hist['bmus'] -1 # para esta funcion bmus tiene que ir de 0 a 41
# data_hist['n_clusters'] = np.arange(42)
#
# data_synth['sorted_bmus'] = data_synth['bmus'] -1 # para esta funcion bmus tiene que ir de 0 a 41
# data_synth['n_clusters'] = np.arange(42)
#
# Plot_DWTs_Mean_Anom_modif(data_hist, data_hist['area'], kind='mean', p_export='/Users/albacid/Projects/SERDP/DWTareas_TWLpeaks_hist.png')
# Plot_DWTs_Mean_Anom_modif(data_synth, data_synth['area'], kind='mean', p_export='/Users/albacid/Projects/SERDP/DWTareas_TWLpeaks_synth.png')
#
# sys.exit()

#----------------------------------------------------------------
# DWT probabilities associated to AWT probabilities

# # Historical
# AWT = ReadMatfile('/Users/albacid/Projects/SERDP/results_files/Historicos/TESLA_AWT.mat')
# #print(AWT['AWT'].keys())
# yy_AWT_hist = np.arange(1880, 2017)
#
# bmus_AWT_hist = []
# for t in data_hist['time']:
#
#     yy_DWT = t.dt.year.values
#     ind = np.where(yy_AWT_hist == yy_DWT)
#     bmus_AWT_hist.append(AWT['AWT']['bmus'][ind][0])
#
#
# bmus_AWT_hist = np.array(bmus_AWT_hist) - 1
# bmus_DWT_hist = data_hist['bmus'] - 1
#
# Plot_Probs_WT_WT(bmus_AWT_hist, bmus_DWT_hist , 6, 42,wt_colors=True,
#                  ttl= 'TWL max DWTs Probabilities by AWTs - Historical', p_export = '/Users/albacid/Projects/SERDP/TWLmax_DWTprob_AWTprob.png')



# Synthetic
AWT = ReadMatfile('/Users/albacid/Projects/SERDP/results_files/Simulados/Hourly_100000years_AWT.mat')

d_ini= np.datetime64('0001-01-01 00:00').astype(datetime)
d_end = np.datetime64('9961-01-06 00:00').astype(datetime)
# d_end = np.datetime64('00012-01-06 00:00').astype(datetime)
time = [d_ini + timedelta(hours=i) for i in range((d_end-d_ini).days * 24)]
yy_AWT_synth = [x.year for x in time]

print(len(AWT['AWT_t']))
print(len(yy_AWT_synth))


bmus_AWT_synth = []
for t in data_synth['time']:

    ind = np.where(yy_AWT_synth == t.dt.year.values)

    bmus_AWT_synth.append(AWT['AWT_t'][ind[0]][-1])# cojo el ultimo valor, porque el primero es diferente durante 1 mes aprox


bmus_AWT_synth = np.array(bmus_AWT_synth) - 1
bmus_DWT_synth = data_synth['bmus'] - 1

print('plotting')

Plot_Probs_WT_WT(bmus_AWT_synth, bmus_DWT_synth , 6, 42,wt_colors=True,
                 ttl= 'TWL max DWTs Probabilities by AWTs - Synthetic', p_export = '/Users/albacid/Projects/SERDP/TWLmax_DWTprob_AWTprob_synth.png')



sys.exit()

#----------------------------------------------------------------
# DWT probabilities associated to MJO fases
MJO_hist = xr.open_dataset('/Users/albacid/Projects/TeslaKit_projects/databases/MJO/MJO_hist.nc')
MJO_phase = MJO_hist.phase.values[:]
MJO_hist['phase'] = MJO_hist.phase - 1

print(data_hist)
print()
print(MJO_hist)
sys.exit()

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


