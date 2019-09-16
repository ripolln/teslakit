import sys
import numpy as np
import xarray as xr
from scipy.io.matlab.mio5_params import mat_struct
import scipy.io as sio
import matplotlib.pyplot as plt

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
#--------------------------------------------------------------------------------
grid = 'GUAM'
x_lim = [144, 145]
y_lim = [13, 14]

lon_GOW = 144.5
lat_GOW = 13.5
lon_CSIRO = 144.6
lat_CSIRO = 13.466

# Load coastfile
coastfile = '/Users/albacid/DATABASES/coastlines/LineaCostaGlobal.mat'
coast = ReadMatfile(coastfile)
lon = coast['costa']['lon']
lat = coast['costa']['lat']

# Load bati GUAM
ds_bati = xr.open_dataset('/Users/albacid/DATABASES/Bathymetries/Guam/depth.nc')


#-----------------------------------------------------------------------------------------------------

# plot bati + GOW + CSIRO

fig = plt.figure(figsize=(12, 9))
c = plt.pcolor(ds_bati['lon'], ds_bati['lat'], ds_bati['elevation'], cmap='gist_ncar', vmin=-4000, vmax=0)
fig.colorbar(c)
plt.plot(lon, lat, 'k')
plt.plot(lon_GOW, lat_GOW, '*r', markersize=10)
plt.plot(lon_CSIRO, lat_CSIRO,'*b', markersize=10)
plt.xlim(x_lim)
plt.ylim(y_lim)
plt.legend(['','GOW','CSIRO'])
plt.title(grid + ' Waves partitions')
plt.gca().set_aspect('equal')
fig.savefig('/Users/albacid/Software/Bitbucket_repos/teslakit/scripts/points_GOW_CSIRO.png', dpi=600)
plt.close()


