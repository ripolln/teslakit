#!/usr/bin/env python
# -*- coding: utf-8 -*-

import xarray as xr
import os.path as op

from lib.objs.predictor import WeatherPredictor as WP
from lib.custom_stats import ClassificationKMA

# data storage
p_data = '/Users/ripollcab/Projects/TESLA-kit/teslakit/data/'


# -------------------------------------------------------------------
# LOAD TESLAKIT PREDICTOR AND DO PRINCIPAL COMPONENTS ANALYSIS

# Load a WeatherPredictor object from netcdf
p_pred_save = op.join(p_data, 'TKPRED_SST.nc')
wpred = WP(p_pred_save)


# calculate running average grouped by months and save
#wpred.CalcRunningMean(5)
#wpred.SaveData()


# Principal Components Analysis
y1 = 1880
yN = 2016
m1 = 6
mN = 5

ipca = wpred.CalcPCA(y1, yN, m1, mN)

#print ipca.components_ (eofs)
#print ipca.explained_variance_

# KMA Classification 
num_clusters = 6
num_reps = 2000
repres = 0.95

AWT = ClassificationKMA(ipca, num_clusters, num_reps, repres)


# ---------------------------------------------------------------------------
# TODO: ENLAZAR CON LA REGRESION LOGISTICA