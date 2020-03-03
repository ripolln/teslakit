
# common
import os
import os.path as op
import sys

# pip
import numpy as np
import xarray as xr
import pandas as pd
np.set_printoptions(precision=2)
np.set_printoptions(suppress=True)

# teslakit
sys.path.insert(0, op.join(op.dirname(__file__),  '..'))
from teslakit.kma import ChangeProbabilities, ClusterProbabilities


rut = '/Users/albacid/Projects/TeslaKit_projects/sites/ROI/SST/'



# historical WT probs
SST_AWTs = xr.open_dataset(rut + 'SST_KMA.nc')
bmus_his = SST_AWTs.bmus.values[:] + 1

num_clusters = 6
set_values = np.arange(num_clusters) + 1
WT_probs_his = ClusterProbabilities(bmus_his, set_values)*100
print(WT_probs_his)



# future WT probs
WT_probs_new = [13.34, 12.07, 19.51, 12.49, 25.19, 17.40] # must sum 100%
print(WT_probs_new)

# markov=1 WT probs
SST_AWTs_sim = xr.open_dataset(rut + 'CC_SST_AWT_sim.nc')
n_sims = len(SST_AWTs_sim.n_sim)

WT_probs_sim_new = np.zeros((n_sims, num_clusters)) * np.nan
for s in range(n_sims):
    bmus_sim = SST_AWTs_sim['evbmus_sims'].isel(n_sim=s).values
    WT_probs_sim_new[s,:] = ClusterProbabilities(bmus_sim, set_values)*100

WT_probs_sim_new_mean = np.mean(WT_probs_sim_new, axis=0)



print(WT_probs_sim_new_mean - WT_probs_his)
print(WT_probs_new - WT_probs_his)



