

#!/usr/bin/env python
# -*- coding: utf-8 -*-

# common
import os
import os.path as op
import json

# pip
import xarray as xr
import numpy as np

# DEV: override installed teslakit
import sys
sys.path.insert(0, op.join(os.path.abspath(''), '..', '..', '..'))

# teslakit
from teslakit.database import Database



# %%
# --------------------------------------
# Teslakit database

p_data = r'/Users/anacrueda/Documents/Proyectos/TESLA/data'
db = Database(p_data)

# set site
db.SetSite('ROI')

# %%
# --------------------------------------
# load data and set parameters

with open('/Users/anacrueda/Documents/Data/STORMs/data.json') as f:
    TCs_STORM = json.load(f)

print(TCs_STORM)


##
# wave point longitude and latitude
pnt_lon = 167.5
pnt_lat = 9.75

# radius for TCs selection (º)
r1 = 14
r2 = 4
