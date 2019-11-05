

import matplotlib.pyplot as plt
import xarray as xr
import sys
from teslakit.io.matlab import ReadMatfile
from teslakit.custom_dateutils import DateConverter_Mat2Py
# import PeakUtils
from scipy.signal import find_peaks


matfile = '/Users/albacid/Projects/SERDP/results_files/Historicos/KWA_historical_parameters_2016_sep.mat'

data = ReadMatfile(matfile)

# print(data.keys())
time = DateConverter_Mat2Py(data['time'])


# 1) Obtain TWL Offshore

def AWL(hs, tp):
   'Returns Atmospheric Water Level'
   # TODO: tp/1.25 ?
   #return 0.043*(hs*1.56*(tp/1.25)**2)**(0.5)
   return 0.043*(hs*1.56*(tp/1.00)**2)**(0.5)

def TWL(awl,ss,at,mmsl):
   'Returns Total Water Level'
   twl =awl+ss+at+mmsl

   return xr.Dataset(
           {
               'TWL': (('time',), twl),
           },
           coords = {'time': time}       )


awl = AWL(data['hs'], data['tp'])
twl = TWL(awl, data['ss'], data['at'], data['mmsl'])

print(twl)
plt.plot(time, twl['TWL'].values)
plt.show()

sys.exit()

out = find_peaks(twl['TWL'].values, threshold=5, distance=None)

print(out)
