import pandas as pd
import numpy as np
import xarray as xr
from datetime import datetime
import sys
np.set_printoptions(precision=2)
np.set_printoptions(suppress=True)
from teslakit.io.aux_nc import StoreBugXdset
from teslakit.kma import ClusterProbabilities, ChangeProbabilities


# inputs
rut = '/Users/albacid/Projects/TeslaKit_projects/sites/ROI/SST/'

bmus_his = xr.open_dataset(rut + 'SST_KMA.nc')
time_his = bmus_his.time.values
bmus_his = bmus_his.bmus.values + 1


# probabilities of each AWT
num_clusters = 6
set_values = np.arange(num_clusters) + 1
WT_probs = ClusterProbabilities(bmus_his, set_values)*100


# # transistion probability matrix (original)
# trans_prob = ChangeProbabilities(bmus_his, set_values)[1] # matriz probs cambio (markov 1)

# transistion probability matrix (modified)
trans_prob = pd.read_csv(rut + 'CC_transition_prob_matrix.csv')
trans_prob = trans_prob.to_numpy()
trans_prob = trans_prob[:, 1:]

print(trans_prob)
print()
print('row sum')
print(np.sum(trans_prob, axis=0))
print(np.sum(np.sum(trans_prob, axis=0)))
print('col sum')
print(np.sum(trans_prob, axis=1))
print(np.sum(np.sum(trans_prob, axis=1)))
print()
print('WT probs')
print(WT_probs)
print('WT probs from transition prob matrix')
print(np.sum(trans_prob, axis=0) / num_clusters*100)
print()


# data for generating bmus
pca_month_ini = 6
y1_sim = 1700
y2_sim = 2700
time_sim = [datetime(y, pca_month_ini,1) for y in range(y1_sim-1, y2_sim+1)]
n_sims = 10000


# bmus generation
bmus_sim = np.zeros((len(time_sim), n_sims), dtype=int) * np.nan

c = 0 # contador
for s in range(n_sims):

    ind_row = bmus_his[:1] -1  # first bmus state. from 0 to 5

    while c < len(time_sim):

        trans_prob_wt = np.cumsum(trans_prob[ind_row, :])

        nrnd = np.random.rand()
        ind_col = np.where(trans_prob_wt > nrnd)[0][0]

        bmus_sim[c, s] = ind_col + 1

        ind_row = ind_col
        c += 1

    c = 0



# obtain AWTs probabilities of synthetic bmus for each sim
WT_probs_sim = np.zeros((n_sims, num_clusters)) * np.nan

for s in range(n_sims):
    WT_probs_sim[s,:] = ClusterProbabilities(bmus_sim[:,s], set_values)*100


# Obtain mean values of AWT probs from all simulations
WT_probs_sim_mean = np.mean(WT_probs_sim, axis=0)
print('WT probs (mean) for synthetic simulations')
print(WT_probs_sim_mean)

print('----------')
# print('WT probs for each synthetic simulation')
# print(WT_probs_sim)
print('----------')
probs_target = [13.34, 12.07, 19.51, 12.49, 25.20, 17.40] # must sum 100%
print(WT_probs_sim_mean-probs_target)

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



