import sys
import  pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from _datetime import datetime, timedelta
from teslakit.util.time_operations import fast_reindex_hourly
from teslakit.io.matlab import ReadMatfile
from teslakit.util.time_operations import DateConverter_Mat2Py
from teslakit.mjo import MJO_Categories
from teslakit.plotting.estela import Plot_DWTs_Probs, Plot_DWTs_Mean_Anom
from teslakit.plotting.wts import Plot_Probs_WT_WT
from teslakit.util.time_operations import xds_common_dates_daily as xcd_daily
from functions_modif import Plot_DWTs_Mean_Anom_modif, Plot_Probs_WT_MJO

rutin = '/Users/albacid/Projects/SERDP/results_files/'


#----------------------------------------------------------------
# 0) load max TWL dataset
data_hist = xr.open_dataset('/Users/albacid/Projects/SERDP/twl_POT_offshore_historical.nc')

data_synth = xr.open_dataset('/Users/albacid/Projects/SERDP/twl_POT_offshore_synthetic.nc')


#----------------------------------------------------------------
# 1) DWT probabilities associated to twl peaks
print(data_hist)
print()
print(data_synth)
print()

Plot_DWTs_Probs(data_hist['bmus'].values, data_hist['time'].values, 42, p_export='/Users/albacid/Projects/SERDP/DWTprob_TWLpeaks_hist.png')
Plot_DWTs_Probs(data_synth['bmus'].values, data_synth['time'].values, 42, p_export='/Users/albacid/Projects/SERDP/DWTprob_TWLpeaks_synth.png')

sys.exit()


# --------------------------------
# 2) plot area values for each DWT

data_hist['sorted_bmus'] = data_hist['bmus'] -1 # para esta funcion bmus tiene que ir de 0 a 41
data_hist['n_clusters'] = np.arange(42)

data_synth['sorted_bmus'] = data_synth['bmus'] -1 # para esta funcion bmus tiene que ir de 0 a 41
data_synth['n_clusters'] = np.arange(42)

Plot_DWTs_Mean_Anom_modif(data_hist, data_hist['area'], kind='mean', p_export='/Users/albacid/Projects/SERDP/DWTareas_TWLpeaks_hist.png')
Plot_DWTs_Mean_Anom_modif(data_synth, data_synth['area'], kind='mean', p_export='/Users/albacid/Projects/SERDP/DWTareas_TWLpeaks_synth.png')

sys.exit()


#----------------------------------------------------------------
# 3) DWT probabilities associated to AWT probabilities

# Historical
AWT = ReadMatfile('/Users/albacid/Projects/SERDP/results_files/Historicos/TESLA_AWT.mat')
yy_AWT_hist = np.arange(1880, 2017) # years.

# find associated AWT for each max
bmus_AWT_hist = []
for t in data_hist['time']:

    yy_DWT = t.dt.year.values
    ind = np.where(yy_AWT_hist == yy_DWT)
    bmus_AWT_hist.append(AWT['AWT']['bmus'][ind][0])


bmus_AWT_hist = np.array(bmus_AWT_hist) - 1
bmus_DWT_hist = data_hist['bmus'] - 1

Plot_Probs_WT_WT(bmus_AWT_hist, bmus_DWT_hist , 6, 42,wt_colors=True,
                 ttl= 'TWL max DWTs Probabilities by AWTs - Historical',show=True)

p_export = '/Users/albacid/Projects/SERDP/TWLmax_DWTprob_AWTprob.png'
sys.exit()


# Synthetic
AWT = ReadMatfile('/Users/albacid/Projects/SERDP/results_files/Simulados/Hourly_100000years_AWT.mat')

d_ini= np.datetime64('0001-01-01 00:00').astype(datetime)
d_end = np.datetime64('9961-01-06 00:00').astype(datetime)
# d_end = np.datetime64('00012-01-06 00:00').astype(datetime)
time = [d_ini + timedelta(hours=i) for i in range((d_end-d_ini).days * 24)]
yy_AWT_synth = [x.year for x in time]


# find associated AWT for each max
bmus_AWT_synth = []
for t in data_synth['time']: # lentisimo!!

    ind = np.where(yy_AWT_synth == t.dt.year.values)

    bmus_AWT_synth.append(AWT['AWT_t'][ind[0]][-1])# cojo el ultimo valor, porque el primero es diferente durante 1 mes aprox


bmus_AWT_synth = np.array(bmus_AWT_synth) - 1
bmus_DWT_synth = data_synth['bmus'] - 1

print('plotting')

Plot_Probs_WT_WT(bmus_AWT_synth, bmus_DWT_synth , 6, 42,wt_colors=True,
                 ttl= 'TWL max DWTs Probabilities by AWTs - Synthetic', p_export = '/Users/albacid/Projects/SERDP/TWLmax_DWTprob_AWTprob_synth.png')

sys.exit()


#----------------------------------------------------------------
# 4) DWT probabilities associated to MJO fases

# Historical
MJO_hist = xr.open_dataset('/Users/albacid/Projects/TeslaKit_projects/databases/MJO/MJO_hist.nc')

categ, d_rmm_categ = MJO_Categories(MJO_hist['rmm1'], MJO_hist['rmm2'], MJO_hist['phase'])
MJO_hist['categ'] = (('time'),(categ - 1))
MJO_hist = MJO_hist.resample(time='1H').pad()


bmus_DWT_hist = data_hist['bmus'] - 1

bmus_MJO_hist = MJO_hist['categ'].sel(time=bmus_DWT_hist['time'].values)


Plot_Probs_WT_MJO(bmus_MJO_hist, bmus_DWT_hist , 25, 42, wt_colors=True,
                 ttl= 'TWL max DWTs Probabilities by MJO category - Historical', show=False, rows=4, cols=7)

plt.savefig('/Users/albacid/Projects/SERDP/TWLmax_DWTprob_MJOprob.png', dpi=600)
plt.close()
sys.exit()


# Synthetic
MJO_synth = xr.open_dataset('/Users/albacid/Projects/TeslaKit_projects/sites/KWAJALEIN/MJO/MJO_sim.nc')

MJO_synth = MJO_synth.isel(n_sim=0)
MJO_synth = MJO_synth.drop(['rmm1', 'rmm2', 'mjo', 'phase'])

MJO_synth = fast_reindex_hourly(MJO_synth)
MJO_synth['categ'] = (('time'),(MJO_synth['evbmus_sims'] - 1))

bmus_DWT_synth = data_synth['bmus'] - 1

bmus_MJO_synth = MJO_synth['categ'].sel(time=bmus_DWT_synth['time'].values)


Plot_Probs_WT_MJO(bmus_MJO_synth, bmus_DWT_synth , 25, 42, wt_colors=True,
                 ttl= 'TWL max DWTs Probabilities by MJO category - Synthetic', show=True, rows=4, cols=7)

plt.savefig('/Users/albacid/Projects/SERDP/TWLmax_DWTprob_MJOprob_synth.png', dpi=600)
plt.close()
