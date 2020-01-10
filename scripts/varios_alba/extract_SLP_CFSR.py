
from teslakit.io.cfs import ReadSLP
import os


p_db = '/Volumes/CLIMUC/DATOS/CFS/prmsl'
lat1 = 60.5
lat2 = -50
lon1 = 115
lon2 = 279
resample = 4 # coge uno de cada 4 datos

rutout = '/Users/albacid/Projects/TeslaKit_projects/sites/GUAM/ESTELA/'
file_name ='SLP_Guam_{0:.1f}_{1:.1f}_{2:.1f}_{3:.1f}.nc'.format(lon1, lon2, lat1, lat2)

p_save = os.path.join(rutout,file_name)
print(p_save)

slp = ReadSLP(p_db, lat1, lat2, lon1, lon2, resample, p_save)
