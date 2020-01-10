

# common
import os.path as op
import sys
# pip
import xarray as xr
import numpy as np
from datetime import date, timedelta, datetime

#teslakit
from teslakit.custom_dateutils import get_years_months_days, datevec2datetime
from teslakit.waves import TWL_AnnualMaxima
from teslakit.io.matlab import ReadMatfile
from teslakit.estela import Plot_DWTs_Probs

# teslakit funtions changed #
#from .waves import TWL - changed to AWL (atmospheric water level)
def AWL(hs, tp):
    'Returns Total Water Level'
    # TODO: tp/1.25 ?
    #return 0.043*(hs*1.56*(tp/1.25)**2)**(0.5)
    return 0.043*(hs*1.56*(tp/1.00)**2)**(0.5)

def TWL(awl,ss,mmsl):
    
    twl =awl+ss+mmsl
    
    return xr.Dataset(
            {
                'TWL': (('time',), twl),
            },
            coords = {'time': twl.time}
        )


#paths
p_historicos = r'/Users/anacrueda/Documents/Proyectos/SERDP/results_files/Historicos/'
p_sim = r'/Users/anacrueda/Documents/Proyectos/SERDP/results_files/Simulados/'

#HISTORICAL
# read mat files
file = op.join(p_historicos, 'KWA_historical_parameters_2016_sep.mat')

db = ReadMatfile(file)
hs = db['hs']
tp = db['tp']
wdir = db['dir']
bmus = db['bmus']
dates = db['dates']
dates_d = datevec2datetime(dates)
at = db['at']
mmsl = db['mmsl']
ss = db['ss']

file_awt = op.join(p_historicos, 'TESLA_AWT_bmus.mat')

db_awt = ReadMatfile(file_awt)
bmus_AWT = db_awt['bmus']
dates_AWT = db_awt['dates']
dates_AWT_d =  datevec2datetime(dates_AWT)

xds_historical = xr.Dataset(
        {
            'bmus' : ('time', bmus),
            'hs'   : ('time', hs),
            'tp'   : ('time', tp),
            'dir'  : ('time', wdir),
            'ss'   : ('time', ss),
            'at'   : ('time', at),
            'mmsl' : ('time', mmsl),
            
        },
        coords = {'time': dates_d}        
)

xds_historical_awt = xr.Dataset(
        {
            'bmus' : ('time', bmus_AWT),
        },
        coords = {'time': dates_AWT_d}
)


# Obtain TWL annual maxima (offshore)
# change TWL por AWL 
xda_AWL = AWL(xds_historical.hs, xds_historical.tp)
xda_TWL = TWL(xda_AWL,xds_historical.ss,xds_historical.mmsl/1000)

#calculate anual maxima
xds_TWL_AMAX = TWL_AnnualMaxima(xda_TWL)


#obtain associated DWT
xds_max_annual = xds_historical.sel(time=xds_TWL_AMAX.time)

#plot wts prob
p_export = '/Users/anacrueda/Documents/Proyectos/TESLA/ROI/results/prob_annualmax_DWT'
bmus = xds_max_annual['bmus'].values[:]  # index to DWT id
bmus_time = xds_max_annual['time'].values[:]
n_clusters = len(np.unique(xds_historical['bmus'].values[:]))
Plot_DWTs_Probs(bmus, bmus_time, n_clusters, p_export)

p_export = '/Users/anacrueda/Documents/Proyectos/TESLA/ROI/results/prob_DWT'
bmus = xds_historical['bmus'].values[:]
bmus_time = xds_historical['time'].values[:]
Plot_DWTs_Probs(bmus,bmus_time,n_clusters,p_export)

print(max(bmus))

##
# SIMULATED
# read mat files
file = op.join(p_sim, 'Hourly_100000years_Hs_Tp_Dir_Level_v2.mat')

db = ReadMatfile(file)
hs = db['HS']
tp = db['TP']
wdir = db['DIR']
bmus = db['BMUS']
#dates = db['dates']
#dates_d = datevec2datetime(dates)
at = db['AT_t']
mmsl = db['MMSL_t']
ss = db['SS']

# Simulation
num_sims = 1
d1_sim = np.datetime64('2020-01-01 00').astype(datetime)
d2_sim = np.datetime64('7020-01-01 00').astype(datetime)



# simulation dates
#extremadamente lento...
dates_sim = [d1_sim + timedelta(hours=i) for i in range((d2_sim-d1_sim).days*24+1)]

#print(dates_sim[-1])

file_awt = op.join(p_sim, 'Hourly_100000years_AWT.mat')

db_awt = ReadMatfile(file_awt)
bmus_AWT = db_awt['AWT_t']
#dates_AWT = db_awt['dates']
#dates_AWT_d = datevec2datetime(dates_AWT)

xds_simulated = xr.Dataset(
    {
        'bmus': ('time', bmus[:len(dates_sim)]),
        'hs': ('time', hs[:len(dates_sim)]),
        'tp': ('time', tp[:len(dates_sim)]),
        'dir': ('time', wdir[:len(dates_sim)]),
        'ss': ('time', ss[:len(dates_sim)]),
        'at': ('time', at[:len(dates_sim)]),
        'mmsl': ('time', mmsl[:len(dates_sim)]),
        'awt' : ('time', bmus_AWT[:len(dates_sim)])

    },
    coords={'time': dates_sim}
)


# Obtain TWL annual maxima (offshore)
# change TWL por AWL
xda_AWL_sim = AWL(xds_simulated.hs, xds_simulated.tp)
xda_TWL_sim = TWL(xda_AWL_sim, xds_simulated.ss, xds_simulated.mmsl / 1000)

# calculate anual maxima
xds_TWL_AMAX_sim = TWL_AnnualMaxima(xda_TWL_sim)

# obtain associated DWT
xds_max_annual_sim = xds_simulated.sel(time=xds_TWL_AMAX_sim.time)

# plot wts prob
# Plot DWTs mean using var_data
p_export_sim = '/Users/anacrueda/Documents/Proyectos/TESLA/ROI/results/prob_annualmax_DWT_simulated'
bmus_sim = xds_max_annual_sim['bmus'].values[:]  # index to DWT id
bmus_time_sim = xds_max_annual_sim['time'].values[:]
n_clusters_sim = len(np.unique(xds_simulated['bmus'].values[:]))
Plot_DWTs_Probs(bmus_sim, bmus_time_sim, n_clusters_sim, p_export_sim)
