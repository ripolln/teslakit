# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 15:48:12 2020

@author: lcag075
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.legend_handler import HandlerBase
import h5py
import xarray as xr
from scipy.stats import circmean
from IPython.display import HTML, display
def print_header(string, string2='', nchar=80):
    print('\n{}\n{}:\n{}\n{}'.format(nchar*'=', string, nchar*'=', string2))

def SuperSpectrum_10degrees(waves, Sectors,st_order):
    cont=0
    for st in st_order: 
        
        pos=((waves.direction.values > Sectors[cont][0]) & (waves.direction.values <= Sectors[cont][1]))

        w1= waves.sel(site=waves.site[st])
        ds = w1.isel(direction = pos)
        ds=ds.drop('site')
        ds=ds.drop('lon')
        ds=ds.drop('lat')
        
        if cont==0:
            ds1=ds;
            del ds
        else:
            ds2=ds1.merge(ds)
            ds1=ds2
            del ds2
            del ds
        
        cont=cont+1
       

    return ds1


def SuperSpectrum_6points(waves, Sectors,st_order):
    cont=0
    for st in st_order: 
        if st == 5:
            pos = (waves.direction.values<30) | (waves.direction.values>330)
        else:
            pos=((waves.direction.values >= Sectors[cont][0]) & (waves.direction.values <= Sectors[cont][1]))

        w1= waves.sel(station=waves.station[st])
        ds = w1.isel(direction = pos)
        ds=ds.drop('station')
        
        if cont==0:
            ds1=ds;
            del ds
        else:
            ds2=ds1.merge(ds)
            ds1=ds2
            del ds2
            del ds
        
        cont=cont+1
       

    return ds1

def PlotSpectrum(dss,time_pos):
    
    fig, ax= plt.subplots(1,figsize=[8,8])
    x=np.deg2rad(dss.direction.values)
    x1=np.append(x,x[0])    
    ds44=dss.isel(time=time_pos,site=0)
    z=ds44.efth.values
    z1=np.column_stack((z[:,:],z[:,-1]))
    z1[np.where(z1<0.000001)]=np.nan
    ax = plt.subplot(111, projection='polar')
    p1=ax.pcolormesh(x1, dss.frequency.values,np.sqrt(z1), vmin=0, vmax=0.2)   
    p1.set_cmap('PuRd')    
    ax.set_theta_zero_location('N', offset=0)
    ax.set_theta_direction(-1)
    cbar = plt.colorbar(p1, ax=ax,pad=0.15, format='%.2f')
    cbar.set_label('Sqrt(Efth)') 
    plt.title('Time: ' + str(time_pos),pad=20)
    
    
    
    
def AddWindEast_Depth(waves_original, p_w):
    

    ds=waves_original.squeeze();
    inif=waves_original.time[0]
    finf=waves_original.time[-1]
    del waves_original
    print(ds)
    unique,unique_ind=np.unique(ds.time,return_index=True) #Hay tiempos repetidos
    ds=ds.isel(time=unique_ind)
    #Redondeamos los minutos 59 a horas enteras    
    ds['time']=ds['time'].dt.round('H').values
        
    # Read wind from file and create xarray dataset
    f = h5py.File(p_w,'r') 
    wspd = f.get('W/speed'); wspd = np.array(wspd) # For converting to numpy array
    wdir = f.get('W/dir') ; wdir = np.array(wdir) # For converting to numpy array
    wdates=f.get('W/dates');    wdates = np.array(wdates) # For converting to numpy array
    depths=np.array([8,4399,4270,3125,3760,3703,3797])
    
    ini=str(int(wdates[0,0]))+'-0'+str(int(wdates[1,0])) +'-0'+str(int(wdates[2,0])) +'T0'+str(int(wdates[3,0])) +':0'+str(int(wdates[4,0]))
    fin=str(int(wdates[0,-1]))+'-0'+str(int(wdates[1,-1])) +'-0'+str(int(wdates[2,-1])) +'T0'+str(int(wdates[3,-1])) +':0'+str(int(wdates[4,-1])) 
    a=np.arange(ini, fin, dtype='datetime64[h]'); wtime=a;
    depths=np.tile(depths, [len(wtime),1])
    
    wspd=wspd[:,0:len(a)]; wdir=wdir[:,0:len(a)];
    stations=[345, 2856, 2857, 2858, 2913, 2914, 2915];
    winds = xr.Dataset({'wdir': (['time','station'],np.transpose(wdir)),'wspd': (['time','station'], np.transpose(wspd)),'dpt': (['time','station'], depths)},coords={'station': stations, 'time': wtime})
    
    winds = winds.sel(time = slice(inif,finf))
    ds['wspd']=winds.wspd.sel(station=2915);  ds['wdir']=winds.wdir.sel(station=2915)
    
    
    depths=[8,4399,4270,3125,3760,3703,3797]
    depth=np.mean(depths[1:]);
    depths=np.tile(depth, [len(ds.time),1]); depths=depths[:,0]
    depths=xr.DataArray(depths,dims=['time'])
    ds['dpt']=depths

    return ds    


def ws(waves, p_w,wcut,msw,agef,model): #WAVESPECTRA: Espectro, pathviento, inicio, fin, estacion, wind cut, max num of swells

       
    unique,unique_ind=np.unique(waves.time,return_index=True) #Hay tiempos repetidos
    waves=waves.isel(time=unique_ind)
    waves['time']=waves['time'].dt.round('H').values#Redondeamos los minutos 59 a horas enteras
    
    # Rename variables 
    
    if model=='spc':
        waves = waves.rename({'frequency':'freq','direction':'dir'}).set_coords({'freq','dir'});
        waves['efth']=waves['efth'] #Los de swan estan en grados
    
    elif model=='csiro':
        waves = waves.rename({'Efth':'efth'})
        waves = waves.rename({'frequency':'freq','direction':'dir'}).set_coords({'freq','dir'});
        waves['efth']=waves['efth']*(np.pi/180); #Para pasar de radianes a grados
        
    ## Watershed partitioning: DEFINITION FUNCTION
    #
    #def partition(self, wsp_darr, wdir_darr, dep_darr, agefac=1.7,
    #                  wscut=0.3333, hs_min=0.001, nearest=False, max_swells=5):
    #   - wsp_darr (DataArray): wind speed (m/s).
    #   - wdir_darr (DataArray): Wind direction (degree).
    #   - dep_darr (DataArray): Water depth (m).
    #   - agefac (float): Age factor.
    #   - wscut (float): Wind speed cutoff.
    #   - hs_min (float): minimum Hs for assigning swell partition.
    #   - nearest (bool): if True, wsp, wdir and dep are allowed to be taken from the.
    #              nearest point if not matching positions in SpecArray (slower).
    #   - max_swells: maximum number of swells to extract
    bulk_params = waves.spec.stats(['hs','tp','tm02','dpm','dspr'])
    ds_part = waves.spec.partition(waves.wspd, waves.wdir, waves.dpt ,wscut=wcut,max_swells=msw, agefac=agef)
    stats_part = ds_part.spec.stats(['hs','tp','tm02','dpm','dspr'])

    return stats_part, ds_part, waves, bulk_params




