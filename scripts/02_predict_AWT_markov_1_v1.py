import pandas as pd
import numpy as np
import xarray as xr
from datetime import datetime
import sys
np.set_printoptions(precision=2)
np.set_printoptions(suppress=True)
from teslakit.io.aux_nc import StoreBugXdset
from teslakit.plotting.wts import Plot_Compare_Transitions
from teslakit.kma import ClusterProbabilities

rut = '/Users/albacid/Projects/TeslaKit_projects/sites/ROI/SST/'

# inputs
trans_prob = pd.read_csv(rut + 'CC_transition_prob_matrix.csv')
trans_prob = trans_prob.to_numpy()
trans_prob = trans_prob[:, 1:]
print(trans_prob)
print(np.shape(trans_prob))

bmus_his = xr.open_dataset(rut + 'SST_KMA.nc')
time_his = bmus_his.time.values
bmus_his = bmus_his.bmus.values


# data for generating bmus
pca_month_ini = 6
y1_sim = 1700
y2_sim = 2700
time_sim = [datetime(y, pca_month_ini,1) for y in range(y1_sim-1, y2_sim+1)]
n_sims = 10


# bmus generation
bmus_sim = np.zeros((len(time_sim), n_sims), dtype=int) * np.nan

c = 0 # contador
for s in range(n_sims):

    ind_row = bmus_his[:1] -1  # first bmus state. from 0 to 5

    while c < len(time_sim):

        trans_prob_wt = np.cumsum(trans_prob[ind_row, :])

        nrnd = np.random.rand()
        ind_col = np.where(trans_prob_wt > nrnd)[0][0]

        #bmus_sim = np.append(bmus_sim, ind_col + 1)
        bmus_sim[c, s] = ind_col + 1

        ind_row = ind_col
        c += 1

    c = 0



# compare TW probabilities with target_probs in 01_modify...
num_clusters = 6
set_values = np.arange(num_clusters) + 1

print('----------')
print(np.sum(trans_prob/6*100, axis=0))
print('----------')
for s in range(n_sims):
    print(ClusterProbabilities(bmus_sim[:,s], set_values)*100)


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



