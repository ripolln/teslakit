
from teslakit.io.cfs import ReadSLP



p_db = '/Volumes/CLIMUC/DATOS/CFS/prmsl'
lat1 = 60.5
lat2 = -50
lon1 = 115
lon2 = 279
resample = 4 # coge uno de cada 4 datos

p_save = '/Users/albacid/Projects/TeslaKit_projects/sites/GUAM/ESTELA'


slp = ReadSLP(p_db, lat1, lat2, lon1, lon2, resample, p_save)
