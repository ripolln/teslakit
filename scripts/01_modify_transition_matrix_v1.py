#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

# input
SST_AWTs = xr.open_dataset(rut + 'SST_KMA.nc')
bmus = SST_AWTs.bmus.values[:]
bmus +=1

num_clusters = 6  # o length permutation, o size m_probs_transition...

probs_target = [13.34, 12.07, 19.51, 12.49, 25.20, 17.40] # must sum 100%
# probs_target = [25.00, 10.00, 30.00, 5.00, 20.00, 10.00]
probs_target = [i/100 for i in probs_target]  # to probabilities instead of percentages

iter_error = 0.0001  # max error in any WT probs


# def custom function to modify probabilities transition matrix by iteration
def permutate_transprobs(trans_prob, target_probs):
    '''
    Permutate probabilities transition matrix, tries to achieve cluster
    probabilities set at probs_target_cluster

    trans_prob - cluster transition probabilities matrix
    target_probs - target permutation cluster probabilities

    returns modified transition probabilities matrix
    '''

    # get probs diff with target
    p_difs = target_probs - np.sum(trans_prob, axis=0) # sum rows

    # modify total probabilities column by column
    for ic, pt in enumerate(p_difs):

        wt_col = trans_prob[:, ic]

        wt_col = wt_col + wt_col/np.sum(wt_col)*pt

        trans_prob[:, ic] = wt_col

    # get probability sums for all rows
    p_sums = np.sum(trans_prob, axis=1)

    # now, modify row by row to make it add to 1/num_clusters
    num_clusters = np.shape(target_probs)[0]
    for ir, ps in enumerate(p_sums):

         trans_prob[ir, :] = trans_prob[ir, :] / (num_clusters*ps)


    return trans_prob


# ---------------------------------------------------------
# before iterate: get target cluster probs (add permutation to base probs)
set_values = np.arange(num_clusters) + 1
m_transition_prob = ChangeProbabilities(bmus, set_values)[1] # matriz probs cambio (markov 1)

m_transition_prob = m_transition_prob / num_clusters

#probs_base = ClusterProbabilities(bmus, set_values)
probs_base = np.sum(m_transition_prob, axis=0)  # sum rows

print(m_transition_prob)
print()
print(probs_base)
print(probs_target)
print(probs_target-probs_base)
print()


# start iterations
x = m_transition_prob.copy()
e = iter_error + 1  # initialize
while e > iter_error:

    # permuate matrix
    x = permutate_transprobs(x, probs_target)

    # check error (max cluster error)
    e = np.max(abs(probs_target - np.sum(x, axis=0)))
    print('error: {0:.5f}'.format(e))

x = x * num_clusters

print('done')
print(x)
print()
print('row sum')
print(np.sum(x, axis=0))
print('col sum')
print(np.sum(x, axis=1))
print()
print(np.sum(np.sum(x, axis=0)))
print(np.sum(np.sum(x, axis=1)))


xds = pd.DataFrame(x)#, columns=['WT1', 'WT2', 'WT3', 'WT4', 'WT5', 'WT6'],
                   #index=['WT1', 'WT2', 'WT3', 'WT4', 'WT5', 'WT6'])

print(xds)
xds.to_csv(rut + 'CC_transition_prob_matrix.csv')
