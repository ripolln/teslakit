

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from datetime import timedelta
import sys
from teslakit.io.matlab import ReadMatfile
from teslakit.custom_dateutils import DateConverter_Mat2Py
from teslakit.waves import TWL as Runup

# read to install scikit-extremes--> https://github.com/kikocorreoso/scikit-extremes
# import skextremes as ske
# No puedo usarlo, POT no está integrado


rutin = '/Users/albacid/Projects/SERDP/results_files/'

#----------------------------------------------------------------
# # 3-hourly offshore values:
# matfile = 'Historicos/KWA_historical_parameters_2016_sep.mat'

# 3-hourly coast values:
matfile = 'Historicos/Reconstr_Hs_Tp_Dir_Level_NOtransectRoi_Historical_2016_p6_in30m_snell.mat'
#----------------------------------------------------------------


# Load data
data = ReadMatfile(matfile)
# print(data.keys())

ss = data['ss']
at = data['at']
mmsl = data['mmsl']/1000.0 # to meters??
hs = data['hs']
tp = data['tp']

time = DateConverter_Mat2Py(data['time'])


#---------------------------------------------------------
# 1) Obtain TWL
runup2 = Runup(hs, tp)

# Offshore: Atmospheric Induced Water level proxy
twl = runup2 + ss  # TWL = AWL

# At the coast: TWL proxy
twl = twl + at + mmsl


fig, ax1 = plt.subplots()
ax1.plot(time, twl, label='twl')


#---------------------------------------------------------
# 2) Obtain threshold value
percentile = 99
twl_threshold = np.percentile (twl, percentile)
print('twl threshold (m): ', round(twl_threshold,2))
print()

ax1.plot([time[0], time[-1]], [twl_threshold, twl_threshold], '-k', label=('threshold:', str(percentile), '%'))


#---------------------------------------------------------
# 3) Obtain duration above threshold and keep maximum peak if consecutive
ind_mask = np.where(twl >= twl_threshold, 1, 0)
ind_dif = np.diff(ind_mask)

ax1.plot(time, twl*ind_mask, '.k', label='twl over threshold')

ind_ini = np.where(ind_dif == 1)[0] + 1
ind_fin = np.where(ind_dif == -1)[0] + 1

time_max = []
twl_max = []
durac = []
for ind_i, ind_f in zip(ind_ini, ind_fin):

    twl_temp = np.max(twl[ind_i:ind_f])
    twl_max.append(twl_temp)

    twl_ind_max = np.argmax(twl[ind_i:ind_f])
    ind_time = list(range(ind_i,ind_f))
    time_temp = time[ind_time[twl_ind_max]]
    time_max.append(time_temp)

    durac.append(time[ind_f]-time[ind_i])

ax1.plot(time_max,twl_max,'*r', label='twl max')


#---------------------------------------------------------
# 4) Ensure independence between maximums
window = timedelta(days=3)

ind_window_mask = np.where(np.diff(time_max) < window, 1, 0)
ind_dif = np.diff(ind_window_mask)

# mover los valores una posicion a la derecha (insertar un 0 en la posición 0)
ind_dif = np.insert(ind_dif, 0, 0)

ind_ini = np.where(ind_dif == 1)[0]
ind_fin = np.where(ind_dif == -1)[0] +1

# Keep maximum of dependent events
time_max_indep = []
twl_max_indep = []
durac_indep = []
ind_delete = []
for ind_i, ind_f in zip(ind_ini, ind_fin):

    twl_temp = np.max(twl_max[ind_i:ind_f])
    twl_max_indep.append(twl_temp)

    twl_ind_max = np.argmax(twl_max[ind_i:ind_f])
    ind_time = list(range(ind_i,ind_f))
    time_temp = time_max[ind_time[twl_ind_max]]
    time_max_indep.append(time_temp)

    durac_indep.append(time_max[ind_f-1]-time_max[ind_i])

    ind_delete.extend(ind_time)


# remove all dependent events
twl_max = np.delete(twl_max, [ind_delete])
time_max = np.delete(time_max, [ind_delete])
durac = np.delete(durac, [ind_delete])

# add independent events
twl_max = np.concatenate((twl_max, twl_max_indep))
time_max = np.concatenate((time_max, time_max_indep))
durac = np.concatenate((durac, durac_indep))
durac = [tt.total_seconds()/3600 for tt in durac] # to hours

# sort data
out = pd.DataFrame({'time': time_max, 'twl': twl_max, 'duration': durac})
out = out.set_index('time')
out.sort_index(inplace=True)

# save
nc = out.to_xarray()
nc.to_netcdf('/Users/albacid/Projects/SERDP/twl_POT.nc')

print('number of extreme events', len(out.twl))
print('max duration of extreme events:', np.max(out.duration[:]), ' h')
n_years = data['dates'][-1][0] - data['dates'][0][0] + 1
print('approximately ', round(len(twl_max)/n_years,1), 'events/yr')


# plot
ax1.plot(out['twl'], '.b', label='twl max indep')

ax2 = ax1.twinx()
ax2.plot(out['duration'], 'og', label='duration')

ax1.set_ylabel('twl (m)',fontdict=dict(weight='bold'))
ax2.set_ylabel('duration (h)',fontdict=dict(weight='bold'))
ax1.legend()
ax2.legend()



plt.show()
