
import os.path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import sys
import scipy.io as sio
from scipy.io.matlab.mio5_params import mat_struct

from datetime import datetime, timedelta
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()


#--------------------------------------------------------------------------------
def matlab_to_python_datetime(datenum_matlab):
    'Return python datetime for matlab datenum. Transform and adjust from matlab.'

    d_new = list(np.ones(len(datenum_matlab),) * np.nan)

    for i, d in enumerate(datenum_matlab):
        d_temp = datetime.fromordinal(int(d)) + \
    timedelta(days=float(d % 1)) - \
        timedelta(days=366) + timedelta(microseconds=0)

        d_new[i] = d_temp

    return d_new

#--------------------------------------------------------------------------------
def ReadMatfile(p_mfile):
    'Parse .mat file to nested python dictionaries'

    def RecursiveMatExplorer(mstruct_data):
        # Recursive function to extrat mat_struct nested contents

        if isinstance(mstruct_data, mat_struct):
            # mstruct_data is a matlab structure object, go deeper
            d_rc = {}
            for fn in mstruct_data._fieldnames:
                d_rc[fn] = RecursiveMatExplorer(getattr(mstruct_data, fn))
            return d_rc

        else:
            # mstruct_data is a numpy.ndarray, return value
            return mstruct_data

    # base matlab data will be in a dict
    mdata = sio.loadmat(p_mfile, squeeze_me=True, struct_as_record=False)
    mdata_keys = [x for x in mdata.keys() if x not in
                  ['__header__','__version__','__globals__']]

    # use recursive function
    dout = {}
    for k in mdata_keys:
        dout[k] = RecursiveMatExplorer(mdata[k])
    return dout


def spatial_gradient(xdset, var_name):
    '''
    Calculate spatial gradient

    xdset:
        (longitude, latitude, time), var_name

    returns xdset with new variable "var_name_gradient"
    '''

    # TODO:check/ ADD ONE ROW/COL EACH SIDE
    var_grad = np.zeros(xdset[var_name].shape)

    Mx = len(xdset.longitude)
    My = len(xdset.latitude)
    lat = xdset.latitude.values

    for it in range(len(xdset.time)):
        var_val = xdset[var_name].isel(time=it).values

        # calculate gradient (matrix)
        m_c = var_val[1:-1,1:-1]
        m_l = np.roll(var_val, -1, axis=1)[1:-1,1:-1]
        m_r = np.roll(var_val, +1, axis=1)[1:-1,1:-1]
        m_u = np.roll(var_val, -1, axis=0)[1:-1,1:-1]
        m_d = np.roll(var_val, +1, axis=0)[1:-1,1:-1]
        m_phi = np.pi*np.abs(lat)/180.0
        m_phi = m_phi[1:-1]

        dpx1 = (m_c - m_l)/np.cos(m_phi[:,None])
        dpx2 = (m_r - m_c)/np.cos(m_phi[:,None])
        dpy1 = m_c - m_d
        dpy2 = m_u - m_c

        vg = (dpx1**2+dpx2**2)/2 + (dpy1**2+dpy2**2)/2
        var_grad[it, 1:-1, 1:-1] = vg

        # calculate gradient (for). old code
        #for i in range(1, Mx-1):
        #    for j in range(1, My-1):
        #        phi = np.pi*np.abs(lat[j])/180.0
        #        dpx1 = (var_val[j,i]   - var_val[j,i-1]) / np.cos(phi)
        #        dpx2 = (var_val[j,i+1] - var_val[j,i])   / np.cos(phi)
        #        dpy1 = (var_val[j,i]   - var_val[j-1,i])
        #        dpy2 = (var_val[j+1,i] - var_val[j,i])
        #        var_grad[it, j, i] = (dpx1**2+dpx2**2)/2 + (dpy1**2+dpy2**2)/2

    # store gradient
    xdset['{0}_gradient'.format(var_name)]= (
        ('time', 'latitude', 'longitude'), var_grad)

    return xdset

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
# data = data.set_index('date')
#
# data = data.replace('NaN', np.nan).astype(float)
#
# ds = data.to_xarray()
# ds.to_netcdf(os.path.join(rutin, 'GUAM_TG_1996_2016.nc'))
#
# sys.exit()


# #------------------------------------------------------------------------
# # 2) Read TG in mat file, save to netcdf
# rutin = '/Users/albacid/Projects/TeslaKit_projects/sites/z_GUAM_Laura/TIDE'
# file = 'Mareografo_Guam.mat'
# TG = ReadMatfile(os.path.join(rutin, file))
#
# time = matlab_to_python_datetime(TG['DATES'])
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

#
# #------------------------------------------------------------------------
# # 3) plot data
# rutin = '/Users/albacid/Projects/TeslaKit_projects'
# TG_sigma = xr.open_dataset(os.path.join(rutin, 'databases', 'Tide', 'Guam','GUAM_TG_1996_2016.nc'))
# print(TG_sigma)

# # plot figure
# fig, axs = plt.subplots(3, 1, figsize=(12,9))
# axs[0].plot(TG_sigma.date, TG_sigma.observed.values, '-r')
# axs[0].plot(TG_sigma.date, TG_sigma.AT, '-b', linewidth = 0.4)
# axs[1].plot(TG_sigma.date, TG_sigma.NTR, '-k', linewidth = 0.4)
# axs[1].plot(TG_sigma.date, TG_sigma.sigma, '-m', linewidth = 0.4)
#
# fig.suptitle('Tide gauge GUAM')
#
# axs[0].legend(['observed SL', 'AT'])
# axs[0].set_xlim([TG_sigma.date[0].values, TG_sigma.date[-1].values])
# axs[0].set_xlabel('time')
# axs[0].set_ylabel('sea level (m)')
# axs[1].legend(['NTR', 'sigma'])
# axs[1].set_xlim([TG_sigma.date[0].values, TG_sigma.date[-1].values])
# axs[1].set_xlabel('time')
# axs[1].set_ylabel('m')


# #------------------------------------------------------------------------
# # Remove outliers in sigma
# TG_sigma['sigma'] = TG_sigma['sigma'].where(TG_sigma['sigma'].values < .8)
#
# d1 = np.datetime64('2011-09-01')
# d2 = np.datetime64('2013-08-31')
# dates = TG_sigma.date.values
#
# p1 = np.where(dates >= d1)[0][0]
# p2 = np.where(dates <= d2)[0][-1]
#
# TG_sigma['sigma'][p1:p2] = np.nan
#
# xds_Sigma = xr.Dataset({ }, coords = {'time': dates},)
# xds_Sigma['sigma'] = (('time',), TG_sigma['sigma'])
# xds_Sigma.to_netcdf(os.path.join('/Users/albacid/Software/Bitbucket_repos/teslakit/scripts','Sigma.nc'))
#
# print(xds_Sigma)


# axs[2].plot(TG_sigma.date, TG_sigma.sigma, '-m', linewidth = 0.4)
# axs[2].legend(['sigma no outliers'])
# axs[2].set_xlim([TG_sigma.date[0].values, TG_sigma.date[-1].values])
# axs[2].set_xlabel('time')
# axs[2].set_ylabel('m')
#
# fig.savefig(os.path.join('/Users/albacid/Software/Bitbucket_repos/teslakit/scripts', 'TG_Guam.png'), dpi=600)


# #------------------------------------------------------------------------
# # Compare sigma to Hs, Tp, Hs**0.5*Tp, mod viento
# Waves = ReadMatfile('/Users/albacid/Projects/TeslaKit_projects/sites/GUAM/WAVES/Particiones_gow2_144.5_13.5.mat')
#
# time = matlab_to_python_datetime(Waves['time'])
#
# xds_Waves = xr.Dataset({ }, coords = {'time': time},)
# xds_Waves['Hs'] = (('time',), Waves['hs'])
# xds_Waves['Tp'] = (('time',), Waves['tp'])
# print(xds_Waves)
# xds_Waves = xds_Waves.sel(time = slice('1996-01-01','2016-12-31'))
# xds_Waves.to_netcdf(os.path.join('/Users/albacid/Software/Bitbucket_repos/teslakit/scripts', 'Waves.nc'))
# sys.exit()

SLP = xr.open_dataset('/Users/albacid/Projects/TeslaKit_projects/databases/ESTELA/Guam/SLP.nc')

# calculate daily gradients
SLP = SLP.sel(time = slice('1996-01-01','2016-12-31'), longitude=slice(150, 180), latitude=slice(20, 0))

xds_SLP_day = spatial_gradient(SLP, 'SLP')

xds_SLP_day = xds_SLP_day.sel(longitude=167, latitude=10.5)
xds_SLP_day = xds_SLP_day.resample(time='1H').mean()

print(xds_SLP_day)
xds_SLP_day.to_netcdf(os.path.join('/Users/albacid/Software/Bitbucket_repos/teslakit/scripts', 'slpGRAD.nc'))



