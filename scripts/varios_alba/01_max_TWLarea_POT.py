

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import xarray as xr
from datetime import timedelta, datetime
import sys
from teslakit.io.matlab import ReadMatfile
from teslakit.custom_dateutils import DateConverter_Mat2Py
from teslakit.waves import AWL as Runup
from teslakit.io.aux_nc import StoreBugXdset


plt.rc('xtick', labelsize=5)
plt.rc('ytick', labelsize=5)

# read to install scikit-extremes--> https://github.com/kikocorreoso/scikit-extremes
# import skextremes as ske
# No puedo usarlo, POT no está integrado

rutin = '/Users/albacid/Projects/SERDP/results_files/'

#----------------------------------------------------------------
# OFFSHORE

# 3-hourly historical offshore values (GOW data):
id = ['offshore', 'historical']
matfile = rutin + 'Historicos/KWA_historical_parameters_2016_sep.mat'

# # # hourly simulated offshore values (output teslakit):
# id = ['offshore', 'synthetic']
# matfile = rutin + 'Simulados/Hourly_100000years_Hs_Tp_Dir_Level_v2.mat'

#----------------
# NEARSHORE

# 3-hourly historical reconstructed values:
#id = ['nearshore', 'historical']
#matfile = rutin + 'Historicos/Reconstr_Hs_Tp_Dir_Level_NOtransectRoi_Historical_2016_p6.mat'

# # hourly simulated reconstructed values:
# id = ['nearshore', 'synthetic']
# matfile = rutin + 'Simulados/Reconstr_Hs_Tp_Dir_Level_NOtransectRoi_100000years_p6_clean.mat'
#----------------------------------------------------------------


#----------------------------------------------------------------
# Load data
data = ReadMatfile(matfile)


if id[0] == 'offshore':
    if id[1] == 'historical':
        ss = data['ss']
        at = data['at']
        mmsl = data['mmsl']/1000.0 # to meters??
        hs = data['hs']
        tp = data['tp']
        bmus = data['bmus']

        time = DateConverter_Mat2Py(data['time'])
        dt_hours = (time[2]-time[1]).total_seconds()/3600.0


    elif id[1] == 'synthetic':

        ss = data['SS']
        at = data['AT_t']
        mmsl = data['MMSL_t']/1000.0 # to meters
        hs = data['HS']
        tp = data['TP']
        bmus = data['BMUS']

        # d_ini= np.datetime64('0001-01-01 00:00')
        # #d_end = np.datetime64('9961-01-05 23:00')
        # time = np.arange(d_ini, d_end + np.timedelta64(1,'h'), dtype='datetime64[h]')
        # dt_hours = (time[2]-time[1])
        # dt_hours = dt_hours / np.timedelta64(1, 'h')

        d_ini= np.datetime64('0001-01-01 00:00').astype(datetime)
        d_end = np.datetime64('9961-01-06 00:00').astype(datetime)
        time = [d_ini + timedelta(hours=i) for i in range((d_end-d_ini).days * 24)]
        dt_hours = (time[2]-time[1]).total_seconds()/3600.0

elif id[0] == 'nearshore':

    # Hs, Tp, Dir, SL(incluye AT, MMSL y SS)
    hs = data['results'][:,0]
    tp = data['results'][:,1]
    dir = data['results'][:,2]
    sl = data['results'][:,3]


    d_ini= np.datetime64('1979-01-01 00:00').astype(datetime)
    d_end = np.datetime64('2016-12-31 00:00').astype(datetime)

    time = [d_ini + timedelta(hours=i) for i in range(0, (d_end-d_ini).days * 24 +1, 3)]
    dt_hours = (time[2]-time[1]).total_seconds()/3600.0

    # waterdepth_NO = data['waterdepth_NO']



#---------------------------------------------------------
# 1) Obtain TWL
runup2 = Runup(hs, tp)

if id[0] == 'offshore':
    # Offshore: Atmospheric Induced Water level proxy
    twl = runup2# + ss  # TWL = AWL

elif id[0] == 'nearshore':
    # At the coast: TWL proxy
    # twl = runup2 + ss + at + mmsl
    twl = runup2 + sl

fig, ax1 = plt.subplots()
ax1.plot(time, twl, label='twl')


#---------------------------------------------------------
# 2) Obtain threshold value
percentile = 99
twl_threshold = np.percentile (twl, percentile)

if id[1] == 'synthetic':
    twl_threshold = 1.38 # threshold del historico offshore

print('twl threshold (m): ', round(twl_threshold,2))
print()

ax1.plot([time[0], time[-1]], [twl_threshold, twl_threshold], '-k', label=('threshold:', str(percentile), '%'))


#---------------------------------------------------------
# 3)From consecutive peaks above threshold: select the maximum, keep duration and
# obtain area above threshold
ind_mask = np.where(twl >= twl_threshold, 1, 0)
ind_dif = np.diff(ind_mask)

ax1.plot(time, twl*ind_mask, '.k', label='twl over threshold')

ind_ini = np.where(ind_dif == 1)[0] + 1
ind_fin = np.where(ind_dif == -1)[0] + 1

time_max = []
twl_max = []
area = []
durac = []
for ind_i, ind_f in zip(ind_ini, ind_fin):

    print('select max from consecutive peaks ', ind_i, '/', ind_fin[-1])

    twl_temp = np.max(twl[ind_i:ind_f])
    twl_max.append(twl_temp)

    area_temp = np.sum(twl[ind_i:ind_f] - twl_threshold) * dt_hours
    area.append(area_temp)

    twl_ind_max = np.argmax(twl[ind_i:ind_f])
    ind_time = list(range(ind_i,ind_f))
    time_temp = time[ind_time[twl_ind_max]]
    time_max.append(time_temp)

    durac.append(time[ind_f]-time[ind_i])

ax1.plot(time_max,twl_max,'*r', label='twl max')


# # plot
# ax2 = ax1.twinx()
# durac = [tt.total_seconds()/3600 for tt in durac] # to hours
# #ax2.plot(time_max, durac, 'og', label='duration')
# #ax2.set_ylabel('duration (h)',fontdict=dict(weight='bold'))
# ax2.plot(time_max, area, 'o', color='grey', label='area')
# ax2.set_ylabel('area (m*h)',fontdict=dict(weight='bold'))
# ax1.set_xlim('1981-02-01','1981-03-01')
# plt.grid()
# plt.show()
# sys.exit()

#---------------------------------------------------------
# 4) Ensure independence between maximums
window = timedelta(days=3)

ind_window_mask = np.where(np.diff(time_max) < window, 1, 0)
ind_dif = np.diff(ind_window_mask)

# mover los valores una posicion a la derecha (insertar un 0 en la posición 0)
ind_dif = np.insert(ind_dif, 0, 0)

# arreglar en caso de que el primer evento ya sea dependiente
if ind_dif[1] == -1:
    ind_dif[0] = 1


ind_ini = np.where(ind_dif == 1)[0]
ind_fin = np.where(ind_dif == -1)[0] +1


# Keep maximum of dependent events
time_max_indep = []
twl_max_indep = []
durac_indep = []
area_indep = []
ind_delete = []
for ind_i, ind_f in zip(ind_ini, ind_fin):

    print('keep independent max ', ind_i, '/', ind_fin[-1])

    twl_temp = np.max(twl_max[ind_i:ind_f])
    twl_max_indep.append(twl_temp)

    twl_ind_max = np.argmax(twl_max[ind_i:ind_f])
    ind_time = list(range(ind_i,ind_f))
    time_temp = time_max[ind_time[twl_ind_max]]
    time_max_indep.append(time_temp)

    durac_indep.append((time_max[ind_f-1]-time_max[ind_i]) + timedelta(hours=dt_hours))
    #durac_indep.append(np.sum(durac[ind_i:ind_f]))

    # print((time_max[ind_f-1]-time_max[ind_i]) + timedelta(hours=dt_hours))
    # print(np.sum(durac[ind_i:ind_f]))
    # print()

    area_indep.append(np.sum(area[ind_i:ind_f]))

    ind_delete.extend(ind_time)


# remove all dependent events
twl_max = np.delete(twl_max, [ind_delete])
time_max = np.delete(time_max, [ind_delete])
durac = np.delete(durac, [ind_delete])
area = np.delete(area, [ind_delete])

# add independent events
twl_max = np.concatenate((twl_max, twl_max_indep))
time_max = np.concatenate((time_max, time_max_indep))
durac = np.concatenate((durac, durac_indep))
durac = [tt.total_seconds()/3600 for tt in durac] # to hours
area = np.concatenate((area, area_indep))


#---------------------------------------------------------
# sort data
out = pd.DataFrame({'time': time_max, 'twl': twl_max, 'twl_exceedances': twl_max-twl_threshold, 'duration': durac, 'area': area})
out = out.set_index('time')

out.sort_index(inplace=True)


#---------------------------------------------------------
# Add the value of the correspondent bmus
temp = pd.DataFrame({'time': time, 'bmus': bmus})
temp = temp.set_index('time')

bmus_twl_max = []
for t in out.index:

    bmus_twl_max.append(temp['bmus'].loc[t])

out['bmus'] = bmus_twl_max


#---------------------------------------------------------
# save
nc = out.to_xarray()

p_ncfile = '/Users/albacid/Projects/SERDP/twl_POT_' + id[0] + '_' + id[1] + '.nc'
StoreBugXdset(nc, p_ncfile)

sys.exit()

#---------------------------------------------------------
print('mean extreme events:', np.mean(out.twl[:]), '(m)')
print('max of extreme events:', np.max(out.twl[:]), '(m)')
print()
print('mean duration of extreme events:', np.mean(out.duration[:]), '(h)')
print('max duration of extreme events:', np.max(out.duration[:]), '(h)')
print()
print('mean area of extreme events:', np.mean(out.area[:]), '(m*h)')
print('max area of extreme events:', np.max(out.area[:]), '(m*h)')
print()
print('number of extreme events:', len(out.twl))
n_years = (time[-1] - time[0]).total_seconds() / (3600*24*365.25)
print('approximately ', round(len(twl_max)/n_years,1), 'events/yr')


#---------------------------------------------------------
# Save fig

# plot
ax1.plot(out['twl'], '.b', label='twl max indep')

ax2 = ax1.twinx()
ax1.set_ylabel('twl (m)',fontdict=dict(weight='bold'))

#ax2.plot(out['duration'], 'og', label='duration')
#ax2.set_ylabel('duration (h)',fontdict=dict(weight='bold'))

ax2.plot(out['area'], 'o', color='grey', label='area')
ax2.set_ylabel('area (m*h)',fontdict=dict(weight='bold'))

ax1.set_xlim('1981-02-01','1981-03-01')

ax1.legend(fontsize=7)
ax2.legend(fontsize=7)

#plt.show()
plt.savefig('/Users/albacid/Projects/SERDP/twl_POT_noSS_' + id[0] + '_' + id[1] + '.png', dpi=600)
