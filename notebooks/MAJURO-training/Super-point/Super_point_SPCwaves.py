# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 14:40:25 2020

@author: lcag075
"""

## DIRECTORY

#import os
#os.chdir (r'C:\Users\lcag075\Dropbox\Culebras-uoa\wavespectra')
import os
os.chdir (r'/Users/laurac/Dropbox/MAJURO/Datos_Killo/UNSWAN_Waves')
import sys
sys.path.append(r'C:\Users\lcag075\Dropbox\Culebras-uoa\wavespectra')
#from lib.SWP import plot_divided_spectrum,PlotPartitions, AddWindEast, SuperSpectrum_5degrees, PlotSpectrum, PlotSpectrumPartitions

from lib.SWP_v2 import SuperSpectrum_10degrees, PlotSpectrum, AddWindEast_Depth, ws

from mpl_toolkits.basemap import Basemap


## LIBRARIES#

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from wavespectra import read_swan #Si se quita no funciona wavespectra
from math import  *
import h5py

mlon=171.18
mlat=7.12

#%% EXAMPLE OF ONE TIME JUST FOR COORDINATES ETC

rs=read_swan('Majuro_bnd_201802')

#%%
mlon=171.18
mlat=7.12

fig=plt.figure(figsize=[16,8])
ax=fig.add_axes([0.1,0.1,0.8,0.8])
m = Basemap(llcrnrlon=170.9,llcrnrlat=6.88,urcrnrlon=171.5,urcrnrlat=7.35, resolution='f',projection='merc',lat_0=7.,lon_0=170., area_thresh=0.000001)
m.drawcoastlines()
m.fillcontinents(color='silver', zorder=5)
m.drawmapboundary(fill_color='lightcyan')
m.drawparallels(np.arange(-60,65,0.1),labels=[1,1,0,1],linewidth=0.2,alpha=0.5)
m.drawmeridians(np.arange(-180,180,0.1),labels=[1,1,0,1],linewidth=0.2, alpha=0.5)
x1,y1=m(mlon,mlat)

NUM_COLORS = int((len(rs.lon)/4))


cm = plt.get_cmap('Paired')
a=[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)];
cm = plt.get_cmap('Set3')
b=[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)];
cm = plt.get_cmap('Accent')
c=[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)];
cm = plt.get_cmap('Set2')
d=[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)];
a=np.concatenate((a,b,c,d))
np.random.shuffle(a) #Changes the order in a

#We calculate the azimuth of the different hindcast points to find the correspondance to the super spectrum

az=np.full([len(rs.lon.values),1],np.nan)

for n in range(len(rs.lon)):
    xx,yy=m(rs.lon.values[n],rs.lat.values[n])
    m.plot(xx,yy, 'bo', markersize=8, color=a[n])
    dLon = np.deg2rad(rs.lon.values[n]) - np.deg2rad(mlon);
    y = sin(dLon) * cos(np.deg2rad(rs.lat.values[n]));
    x = cos(np.deg2rad(mlat))*sin(np.deg2rad(rs.lat.values[n])) -  sin(np.deg2rad(mlat))*cos(np.deg2rad(rs.lat.values[n]))*cos(dLon);
    brng = np.rad2deg(atan2(y, x));
    if brng < 0: brng+= 360
    az[n]=brng


#We create a vector stating which of the different spectral sites contains the information for each 10 degrees sector

sect_dir=np.arange(5,358,10)
xx,yy=m(rs.lon.values,rs.lat.values)
site=np.full([len(sect_dir),1],np.nan)


for mm in range(len(sect_dir)):    
    d=abs(az-sect_dir[mm])
    site[mm]=np.where(d==np.nanmin(d))[0]
    m.plot(xx[site[mm].astype('int')],yy[site[mm].astype('int')] ,markersize=8, color=a[site[mm][0].astype('int')])
    z=[x1,x1+99999999999999*np.sin(np.deg2rad(sect_dir[mm]))]
    z1=[y1,y1+99999999999999*np.cos(np.deg2rad(sect_dir[mm]))]
    m.plot(z,z1,linewidth=0.8,color=a[site[mm][0].astype('int')])
    
#fig.savefig(r'C:\Users\lcag075\Dropbox\Culebras-uoa\Figures\Super-Point-10degrees.png',dpi=300)

#%% SUPER POINT

st_order=site.astype('int') #Stations and order 2914 2915 2858 2857 2856 2913
a=np.arange(0,351,10)
b=np.arange(10,361,10)
Sectors=np.transpose([a.T,b.T])

rs = rs.rename({'freq':'frequency','dir':'direction'}).set_coords({'frequency','direction'})
#waves_original = xr.Dataset({'Efth': (['time','site','frequency','direction'],rs.efth.values)}, coords={'time': rs.time.values, 'site': rs.site.values, 'frequency': rs.freq.values, 'direction': rs.dir.values})

waves=SuperSpectrum_10degrees(rs, Sectors,st_order)

time_pos=56
PlotSpectrum(waves,time_pos) #Full spectrum


fig, ax= plt.subplots(1,figsize=[8,8])
x=np.deg2rad(rs.direction.values)
x1=np.append(x,x[0])    
ds44=rs.isel(time=time_pos,site=19)
z=ds44.efth.values
z1=np.column_stack((z[:,:],z[:,-1]))
z1[np.where(z1<0.000001)]=np.nan
ax = plt.subplot(111, projection='polar')
p1=ax.pcolormesh(x1, rs.frequency.values,np.sqrt(z1), vmin=0, vmax=0.2)   
p1.set_cmap('PuRd')    
ax.set_theta_zero_location('N', offset=0)
ax.set_theta_direction(-1)
cbar = plt.colorbar(p1, ax=ax,pad=0.15, format='%.2f')
cbar.set_label('Sqrt(Efth)') 
plt.title('Time: ' + str(time_pos),pad=20)

#%%

st_order=site.astype('int') #Stations and order
a=np.arange(0,351,10); b=np.arange(10,361,10); Sectors=np.transpose([a.T,b.T])
#p_wk=r'C:\Users\lcag075\Dropbox\MAJURO\Datos_Killo\ERA5_winds\ERA5_Winds_Majuro.mat'
p_w = os.path.join(r'C:\Users\lcag075\Dropbox\Culebras-uoa\MAJURO/DATA/WIND/', "Winds_Majuro_stations_1980_2018.mat")
#Wavespectra parameters
wcut=0.00000000001 #!!!! Wind cut 
msw=15 #!!!! Max number of swells 
agef=1.7 #!!!! Age Factor
model='spc' #Choose between spc or csiro

Yini=1979
Yend=2019
a=np.arange(Yini,Yend,1)
for ty in a:
    for m in range(12):
        
        rs=read_swan('Majuro_bnd_'+str(ty)+str(m+1).zfill(2))
        rs = rs.rename({'freq':'frequency','dir':'direction'}).set_coords({'frequency','direction'})
        waves=SuperSpectrum_10degrees(rs, Sectors,st_order)
        
#        time_pos=100
#        PlotSpectrum(waves,time_pos) #Full spectrum
        
        waves_wind=AddWindEast_Depth(waves, p_w)
        spec_part, ds_part, ds, bulk_params = ws(waves_wind,p_w,wcut,msw,agef,model)
        ds_part.to_netcdf(path=r'C:\Users\lcag075\Dropbox\Culebras-uoa\Majuro_partitions_SPC/DS_PART_'+ str(ty) + '_' + str(m+1).zfill(2) +'.nc')
        ds.to_netcdf(path=r'C:\Users\lcag075\Dropbox\Culebras-uoa\Majuro_partitions_SPC/DS_'+ str(ty) + '_' + str(m+1).zfill(2) +'.nc')
        
        if m==0:
            PARTITIONS=spec_part
            BULK=bulk_params
        else:
            PARTITIONS=PARTITIONS.merge(spec_part)
            BULK=BULK.merge(bulk_params)
        
    PARTITIONS.to_netcdf(path=r'C:\Users\lcag075\Dropbox\Culebras-uoa\Majuro_partitions_SPC/PARTITIONS_'+ str(ty) +'.nc')    
    BULK.to_netcdf(path=r'C:\Users\lcag075\Dropbox\Culebras-uoa\Majuro_partitions_SPC/BULK_'+ str(ty) +'.nc')


#%% WE JOIN THE DIFFERENT YEARS
    
partitions_dir=r'C:\Users\lcag075\Dropbox\Culebras-uoa\Majuro_partitions_SPC'
 
ds = xr.open_dataset(os.path.join(partitions_dir,'PARTITIONS_' + str(1979) + '.nc'))
ds_b = xr.open_dataset(os.path.join(partitions_dir,'BULK_' + str(1979) + '.nc'))

Yini=1980
Yend=2019
a=np.arange(Yini,Yend,1)    
for ty in a:
    dsp =xr.open_dataset(os.path.join(partitions_dir,'PARTITIONS_' + str(ty) + '.nc'))
    ds=ds.merge(dsp)
    
    dspb =xr.open_dataset(os.path.join(partitions_dir,'BULK_' + str(ty) + '.nc'))
    ds_b=ds_b.merge(dspb)
    
    
ds.to_netcdf(path=r'C:\Users\lcag075\Dropbox\Culebras-uoa\Majuro_partitions_SPC/PARTITIONS_1979_2018.nc')
ds_b.to_netcdf(path=r'C:\Users\lcag075\Dropbox\Culebras-uoa\Majuro_partitions_SPC/BULK_1979_2018.nc')
      
