
from datetime import datetime
import numpy as np
np.set_printoptions(precision=2)
np.set_printoptions(suppress=True)
import xarray as xr
import sys
from teslakit.kma import ClusterProbabilities
from teslakit.io.aux_nc import StoreBugXdset


rut = '/Users/albacid/Projects/TeslaKit_projects/sites/ROI/SST/'


# Load historical AWT time series:
SST_AWTs = xr.open_dataset(rut + 'SST_KMA.nc')
bmus = SST_AWTs.bmus.values[:] + 1


# Obtain probability of each AWT
num_clusters = 6
set_values = np.arange(num_clusters) + 1
WT_probs = ClusterProbabilities(bmus, set_values)*100
print(WT_probs)


# Obtain expected change
expected_dif = [2.20, -7.56, 3.49, -9.38, 24.10, -14.29] # From Mathew
WT_probs_new = WT_probs + WT_probs*expected_dif/100
WT_probs_new = WT_probs_new/np.sum(WT_probs_new)
print(WT_probs_new*100)


# Create future bmus with expected probs

# data for generating bmus
pca_month_ini = 6
y1_sim = 1700
y2_sim = 2700
time_sim = [datetime(y, pca_month_ini,1) for y in range(y1_sim-1, y2_sim+1)]
n_sims = 100

# bmus generation
WT_probs_cumsum = np.cumsum(WT_probs_new)
bmus_sim = np.zeros((len(time_sim), n_sims), dtype=int) * np.nan

for s in range(n_sims):

    c = 0
    while c < len(time_sim):

        nrnd = np.random.rand()
        ind_col = np.where(WT_probs_cumsum > nrnd)[0][0]

        bmus_sim[c, s] = ind_col + 1
        c += 1


print('----------')
WT_probs_sim = np.zeros((n_sims, num_clusters)) * np.nan

for s in range(n_sims):
    WT_probs_sim[s,:] = ClusterProbabilities(bmus_sim[:,s], set_values)*100

print(WT_probs_sim)
print('----------')

WT_probs_sim_mean = np.mean(WT_probs_sim, axis=0)
print(WT_probs_sim_mean)

sys.exit()

# Save
xds_bmus_sim = xr.Dataset(
    {
        'evbmus_sims' : (('time','n_sim'), bmus_sim.astype(int)),
    },
    coords = {
        'time' : time_sim,
        'n_sim' : range(n_sims),
    },
)

print(xds_bmus_sim)


StoreBugXdset(xds_bmus_sim, rut + 'CC_SST_AWT_sim.nc')


