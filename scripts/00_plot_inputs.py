

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



rutin = '/Users/albacid/Projects/SERDP/results_files/'

#----------------------------------------------------------------
# OFFSHORE

# # 3-hourly historical offshore values (GOW data):
# id = ['offshore', 'historical']
# matfile = rutin + 'Historicos/KWA_historical_parameters_2016_sep.mat'

# hourly simulated offshore values (output teslakit):
id = ['offshore', 'synthetic']
matfile = rutin + 'Simulados/Hourly_100000years_Hs_Tp_Dir_Level_v2.mat'

#----------------
# NEARSHORE

# # 3-hourly historical reconstructed values:
# id = ['nearshore', 'historical']
# matfile = rutin + 'Historicos/Reconstr_Hs_Tp_Dir_Level_NOtransectRoi_Historical_2016_p6.mat'


# # hourly simulated reconstructed values:
# id = ['nearshore', 'synthetic']
# matfile = rutin + 'Simulados/Reconstr_Hs_Tp_Dir_Level_NOtransectRoi_100000years_p6_clean.mat'
#----------------------------------------------------------------


#----------------------------------------------------------------
# Load data
data = ReadMatfile(matfile)

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

    time = np.arange(0, len(hs))

# # plot components
# fig, ax = plt.subplots(4)
#
# ax[0].plot(time, ss)
# ax[0].set_ylabel('ss (m)')
# ax[1].plot(time, at)
# ax[1].set_ylabel('at (m)')
# ax[2].plot(time, ss)
# ax[2].set_ylabel('mmsl (m)')
# ax[3].plot(time, ss)
# ax[3].set_ylabel('hs (m)')
#
# plt.savefig('/Users/albacid/Projects/SERDP/inputdata_' + id[0] + '_' + id[1] + '.png', dpi=600)


# plot TWL
runup2 = Runup(hs, tp)

if id[0] == 'offshore':
    # Offshore: Atmospheric Induced Water level proxy
    twl = runup2 + ss  # TWL = AWL

elif id[0] == 'nearshore':
    # At the coast: TWL proxy
    twl = runup2 + ss + at + mmsl


fig, ax1 = plt.subplots()
ax1.plot(time, twl)
ax1.plot([time[0], time[-1]], [0.97, 0.97])
ax1.set_ylabel('twl ' + id[0])
plt.show()
sys.exit()
plt.savefig('/Users/albacid/Projects/SERDP/inputdata_twl_' + id[0] + '_' + id[1] + '.png', dpi=300)
