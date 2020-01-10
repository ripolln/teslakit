import xarray as xr
import os.path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import pearsonr
import sys


rutin = '/Users/albacid/Software/Bitbucket_repos/teslakit/scripts'

sigma = xr.open_dataset(os.path.join(rutin, 'Sigma.nc'))
waves = xr.open_dataset(os.path.join(rutin, 'Waves.nc')) # 3-h values
slpGrad = xr.open_dataset(os.path.join(rutin, 'slpGRAD.nc'))

# select coincident values
sigma = sigma.sel(time=waves.time.values)
slpGrad = slpGrad.sel(time=waves.time.values)

time = waves.time.values

# normalize variables
sigma = (sigma['sigma'].values - np.nanmin(sigma['sigma'].values)) / (np.nanmax(sigma['sigma'].values) - np.nanmin(sigma['sigma'].values))
hs = (waves['Hs'].values - np.nanmin(waves['Hs'].values)) / ( np.nanmax(waves['Hs'].values) -  np.nanmin(waves['Hs'].values))
tp = (waves['Tp'].values - np.nanmin(waves['Tp'].values)) / ( np.nanmax(waves['Tp'].values) -  np.nanmin(waves['Tp'].values))
grad = (slpGrad['SLP_gradient'].values - np.nanmin(slpGrad['SLP_gradient'].values)) / (np.nanmax(slpGrad['SLP_gradient'].values) - np.nanmin(slpGrad['SLP_gradient'].values))

# get not nan values
ind = np.where(~np.isnan(sigma))
time = time[ind]
sigma = sigma[ind]
hs = hs[ind]
tp = tp[ind]
grad = grad[ind]


# figure
fig = plt.figure(figsize=(12,9))
gs = GridSpec(4,4)

# sigma vs hs
var = hs
corr = pearsonr(sigma, var)
pos = 0
ax1 = fig.add_subplot(gs[pos,:3])
ax1.plot(time, sigma, 'r')
ax1.plot(time, var)
ax1.legend(['sigma', 'hs'])
ax2 = fig.add_subplot(gs[pos,3])
ax2.plot(sigma, var, '.', markersize=0.5)
ax2.plot([0, 1],[0, 1], color=[.5,.5,.5])
ax2.set_aspect('equal')
ax2.set_ylabel('hs')
ax2.text(.1,.9, 'corr={0:.1f}'.format(corr[0]))


# sigma vs tp
var = tp
corr = pearsonr(sigma, var)
pos = 1
ax1 = fig.add_subplot(gs[pos,:3])
ax1.plot(time, sigma, 'r')
ax1.plot(time, var, 'b')
ax1.legend(['sigma', 'tp'])
ax2 = fig.add_subplot(gs[pos,3])
ax2.plot(sigma, var, '.', markersize=0.5)
ax2.plot([0, 1],[0, 1], color=[.5,.5,.5])
ax2.set_aspect('equal')
ax2.set_ylabel('tp')
ax2.text(.1,.9, 'corr={0:.1f}'.format(corr[0]))

# sigma vs hs^.5*Tp (normalized)
var = (waves['Hs'].values**0.5)*waves['Tp'].values
var = (var - np.nanmin(var)) / (np.nanmax(var) - np.nanmin(var))
var = var[ind]
corr = pearsonr(sigma, var)
pos = 2
ax1 = fig.add_subplot(gs[pos,:3])
ax1.plot(time, sigma, 'r')
ax1.plot(time, var, 'g')
ax1.legend(['sigma', 'hs^.5*tp'])
ax2 = fig.add_subplot(gs[pos,3])
ax2.plot(sigma, var, '.', markersize=0.5)
ax2.plot([0, 1],[0, 1], color=[.5,.5,.5])
ax2.set_aspect('equal')
ax2.set_ylabel('Hs^.5*Tp')
ax2.text(.1,.9, 'corr={0:.1f}'.format(corr[0]))

# sigma vs wind (grad SLP)
var = grad
corr = pearsonr(sigma, var)
pos = 3
ax1 = fig.add_subplot(gs[pos,:3])
ax1.plot(time, sigma, 'r')
ax1.plot(time, var, 'y')
ax1.legend(['sigma', 'gradSLP'])
ax2 = fig.add_subplot(gs[pos,3])
ax2.plot(sigma, var, '.', markersize=0.5)
ax2.plot([0, 1],[0, 1], color=[.5,.5,.5])
ax2.set_aspect('equal')
ax2.set_ylabel('gradSLP')
ax2.set_xlabel('sigma')
ax2.text(.1,.9, 'corr={0:.1f}'.format(corr[0]))

fig.suptitle('Normalized variables GUAM')

fig.savefig(os.path.join('/Users/albacid/Software/Bitbucket_repos/teslakit/scripts', 'Sigma_Guam.png'), dpi=600)


