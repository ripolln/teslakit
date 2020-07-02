#!/usr/bin/env python
# -*- coding: utf-8 -*-

from math import radians, degrees, sin, cos, asin, acos, sqrt, atan2, pi
import numpy as np
import xarray as xr
import os.path as op
import glob
import pandas as pd

from .util.operations import GetUniqueRows


def ReadStormFiles(p_storm, basins):
    '''
    Read STORM files for given basin
    Return xarray.Dataset
    '''

    ls_path = op.join(p_storm, 'STORM_DATA_IBTRACS_{}*.txt'.format([basin for basin in basins]))

    file_list = sorted(glob.glob(ls_path))

    tc_list = []
    cont = 0
    for file in file_list:
        df = pd.read_table(file,
                           sep=',',
                           header=None,
                           names=(
                           'year', 'month', 'storm', 'TimeStep', 'Basin', 'lat', 'lon', 'pressure_min', 'windspeed_max',
                           'rmax', 'cat', 'land', 'dist_land'),
                           dtype={'year': np.int64,
                                  'month': np.int64,
                                  'storm': np.int64,
                                  'TimeStep': np.int64,
                                  'basin': np.int64,
                                  'lat': np.float64,
                                  'lon': np.float64,
                                  'pressure_min': np.float64,
                                  'windspeed_max': np.float64,
                                  'rmax': np.float64,
                                  'cat': np.int64,
                                  'land': np.int64,
                                  'dist_land': np.float64})

        for i in df['year'].unique():
            if i == 0: continue
            late_storm = np.nanmax(df.loc[df.year == i - 1, ['storm']])
            df.loc[df.year == i, ['storm']] = df.loc[df.year == i, ['storm']] + late_storm + 1

        if cont > 0:
            last_storm = storms[-1]
            last_year = years[-1]
            df['storm'] = df['storm'] + last_storm + 1
            df['year'] = df['year'] + last_year + 1
        df.set_index('storm', inplace=True)
        # store DataFrame in list
        tc_list.append(df)
        storms = df.index.tolist()
        years = df.year.tolist()
        cont = cont + 1

    tc_list = pd.concat(tc_list)
    tc_storms = tc_list.to_xarray()

    return tc_storms


def GeoDistance(lat1, lon1, lat2, lon2):
    'Returns great circle distance between points in degrees'

    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    a = sin((lat2-lat1)/2)**2 + cos(lat1) * cos(lat2) * sin((lon2-lon1)/2)**2;
    if a < 0: a = 0
    if a > 1: a = 1

    r = 1
    rng = r * 2 * atan2(sqrt(a), sqrt(1-a))
    rng = degrees(rng)

    return rng

def GeoAzimuth(lat1, lon1, lat2, lon2):
    'Returns geodesic azimuth between point1 and point2'

    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    az = atan2(
        cos(lat2) * sin(lon2-lon1),
        cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(lon2-lon1)
    )
    if lat1 <= -pi/2: az = 0
    if lat2 >=  pi/2: az = 0
    if lat2 <= -pi/2: az = pi
    if lat1 >=  pi/2: az = pi

    az = az % (2*pi)
    az = degrees(az)

    return az
    
def Extract_Circle_STORM(xds_TCs, p_lon, p_lat, r, d_vns):
    '''
    Extracts TCs inside circle - used with STORM database

    xds_TCs: tropical cyclones track database
        lon, lat, pressure variables
        storm dimension

    circle defined by:
        p_lon, p_lat  -  circle center
        r             -  circle radius (degree)

    d_vns: dictionary to set longitude, latitude, time and pressure varnames

    returns:
        xds_area: selection of xds_TCs inside circle
        xds_inside: contains TCs custom data inside circle
    '''

    # point longitude and latitude
    lonlat_p = np.array([[p_lon, p_lat]])

    # get names of vars: longitude, latitude, pressure and time
    nm_lon = d_vns['longitude']
    nm_lat = d_vns['latitude']
    nm_prs = d_vns['pressure']
    nm_rad = d_vns['radius']
    nm_win = d_vns['mwinds']
    nm_tim = d_vns['time']

    # storms longitude, latitude, pressure and time (if available)
    lon = xds_TCs[nm_lon].values[:]
    lat = xds_TCs[nm_lat].values[:]
    prs = xds_TCs[nm_prs].values[:]
    radpm = xds_TCs[nm_rad].values[:]
    mwin = xds_TCs[nm_win].values[:]
    time = xds_TCs[nm_tim].values[:]
    storms =xds_TCs['storm'].values[:]
    # get storms inside circle area
    n_storms = len(np.unique(xds_TCs.storm))
    max_time = np.max(time)
    #print('num storms')
    #print(n_storms)

    #print('storms')
    #print(storms)

    l_storms_area = []

    # inside parameters holders
    l_prs_min_in = []   # circle minimun pressure
    l_prs_mean_in = []  # circle mean pressure
    l_vel_mean_in = []  # circle mean translation speed
    l_rad_mean_in = []  # circle mean radius max winds
    l_wind_mean_in = []
    l_categ_in = []     # circle storm category
    l_date_in = []      # circle date (day)
    l_date_last = []    # last cyclone date 
    l_gamma = []        # azimuth 
    l_delta = []        # delta 

    l_ix_in = []        # historical enters the circle index
    l_ix_out = []       # historical leaves the circle index

    print(np.shape(np.meshgrid(np.ones(max_time) * np.nan, np.ones(n_storms) * np.nan)[0]))
    # dataset para almacenar las tracks
    xds_TCs_sel = xr.Dataset(
        {
            'pressure': (('storm', 'time'), np.meshgrid(np.ones(max_time) * np.nan, np.ones(n_storms) * np.nan)[0]),
            'velocity': (('storm', 'time'), np.meshgrid(np.ones(max_time) * np.nan, np.ones(n_storms) * np.nan)[0]),
            'wind': (('storm', 'time'), np.meshgrid(np.ones(max_time) * np.nan, np.ones(n_storms) * np.nan)[0]),
            'radius': (('storm', 'time'), np.meshgrid(np.ones(max_time) * np.nan, np.ones(n_storms) * np.nan)[0]),
            'lon': (('storm', 'time'), np.meshgrid(np.ones(max_time) * np.nan, np.ones(n_storms) * np.nan)[0]),
            'lat': (('storm', 'time'), np.meshgrid(np.ones(max_time) * np.nan, np.ones(n_storms) * np.nan)[0]),
        },
        coords={
            'storm': (np.arange(n_storms)),
            'time': (np.arange(max_time))
        }
    )
    cont = 0
    for i_storm in range(n_storms):
        # fix longitude <0 data and skip "one point" tracks
        # pos_storm = [i for i, value in enumerate(storms) if value == i_storm]
        pos_storm = [storms == i_storm]
        # print('tormenta i')
        # print(i_storm)
        # print('posicion_tormenta i')
        # print(pos_storm)
        lon_storm = lon[pos_storm]
        #print(lon_storm)
        if not isinstance(lon_storm, np.ndarray): continue
        lon_storm[lon_storm<0] = lon_storm[lon_storm<0] + 360
        #print('longitude')
        #print(lon_storm)

        # stack storm longitude, latitude
        lonlat_s = np.column_stack(
            (lon_storm, lat[pos_storm])
        )

        # index for removing nans
        ix_nonan = ~np.isnan(lonlat_s).any(axis=1)
        lonlat_s = lonlat_s[ix_nonan]

        # calculate geodesic distance (degree)
        geo_dist = []
        for lon_ps, lat_ps in lonlat_s:
            geo_dist.append(GeoDistance(lat_ps, lon_ps, p_lat, p_lon))
        geo_dist = np.asarray(geo_dist)

        # find storm inside circle and calculate parameters
        if (geo_dist < r).any() & (geo_dist.size > 1):
            # print(cont)
            # storm inside circle
            ix_in = np.where(geo_dist < r)[0][:]

            # storm translation velocity
            geo_dist_ss = []
            for i_row in range(lonlat_s.shape[0]-1):
                i0_lat, i0_lon = lonlat_s[i_row][1], lonlat_s[i_row][0]
                i1_lat, i1_lon = lonlat_s[i_row+1][1], lonlat_s[i_row+1][0]
                geo_dist_ss.append(GeoDistance(i0_lat, i0_lon, i1_lat, i1_lon))
            geo_dist_ss = np.asarray(geo_dist_ss)

            # print('pos_storm')
            # print(pos_storm)

            # get delta time in hours (STORM database 3hourly time-steps)
            delta_h = 3 #len(time[pos_storm]) #*3

            #print('time')
            #print(time[pos_storm])

            #print('delta_h')
            #print(delta_h)

            #print('geo_dist_ss')
            #print(geo_dist_ss * 111.0)

            vel = geo_dist_ss * 111.0/delta_h  # km/h

            #print('vel')
            #print(vel)

            # promediate vel 
            velpm = (vel[:-1] + vel[1:])/2
            velpm = np.append(vel[0], velpm)
            velpm = np.append(velpm, vel[-1])

            # calculate azimuth 
            lat_in_end, lon_in_end = lonlat_s[ix_in[-1]][1], lonlat_s[ix_in[-1]][0]
            lat_in_ini, lon_in_ini = lonlat_s[ix_in[0]][1], lonlat_s[ix_in[0]][0]
            gamma = GeoAzimuth(lat_in_end, lon_in_end, lat_in_ini, lon_in_ini)
            if gamma < 0.0: gamma += 360

            # calculate delta
            nd = 1000
            st = 2*np.pi/nd
            ang = np.arange(0, 2*np.pi + st, st)
            xps = r * np.cos(ang) + p_lat
            yps = r * np.sin(ang) + p_lon
            angle_radius = []
            for x, y in zip(xps, yps):
                angle_radius.append(GeoAzimuth(lat_in_end, lon_in_end, x, y))
            angle_radius = np.asarray(angle_radius)

            im = np.argmin(np.absolute(angle_radius - gamma))
            delta = GeoAzimuth(p_lat, p_lon, xps[im], yps[im]) # (-180, +180)
            if delta < 0.0: delta += 360

            #print('pressure')
            #print(prs)

            # more parameters 
            prs_s_in = prs[pos_storm][ix_in]  # pressure
            prs_s_min = np.min(prs_s_in)  # pressure minimun
            prs_s_mean = np.mean(prs_s_in)

            vel_s_in = velpm[ix_in]  # velocity
            vel_s_mean = np.mean(vel_s_in)# velocity mean

            rad_s_in = radpm[pos_storm][ix_in] # radius
            rad_s_mean = np.mean(rad_s_in) # mean radius

            wind_s_in = mwin[pos_storm][ix_in] # max winds
            wind_s_mean = np.mean(wind_s_in) # mean max winds

            categ = GetStormCategory(prs_s_min)  # category

            dist_in = geo_dist[ix_in]
            p_dm = np.where((dist_in==np.min(dist_in)))[0]  # closest to point

            time_s_in = time[pos_storm][ix_in]  # time
            time_closest = time_s_in[p_dm][0]  # time closest to point 

            # filter storms 
            # TODO: storms with only one track point inside radius. solve?
            if np.isnan(np.array(prs_s_in)).any() or \
               (np.array(prs_s_in) <= 860).any() or \
               gamma == 0.0:
                continue

            # store parameters
            l_storms_area.append(cont)
            l_prs_min_in.append(np.array(prs_s_min))
            l_prs_mean_in.append(np.array(prs_s_mean))
            l_vel_mean_in.append(np.array(vel_s_mean))
            l_rad_mean_in.append(np.array(rad_s_mean))
            l_wind_mean_in.append((np.array(wind_s_mean)))
            l_categ_in.append(np.array(categ))
            l_date_in.append(time_closest)
            l_gamma.append(gamma)
            l_delta.append(delta)

            # store storm parameters

            xds_TCs_sel.pressure[cont, :len(time[pos_storm])] = prs[pos_storm]
            xds_TCs_sel.velocity[cont, :len(time[pos_storm])] = velpm
            xds_TCs_sel.wind[cont, :len(time[pos_storm])] = mwin[pos_storm]
            xds_TCs_sel.radius[cont,:len(time[pos_storm])] = radpm[pos_storm]
            xds_TCs_sel.lon[cont, :len(time[pos_storm])] = lon_storm
            xds_TCs_sel.lat[cont,:len(time[pos_storm])] = lat[pos_storm]

            # store historical indexes inside circle
            l_ix_in.append(ix_in[0])
            l_ix_out.append(ix_in[-1])
            cont = cont + 1
            # store last cyclone date too
            #l_date_last.append(time[i_storm][ix_nonan][-1])

    # cut storm dataset to selection
    # xds_TCs_sel = xds_TCs.isel(storm = l_storms_area)
    # xds_TCs_sel = xds_TCs_sel.assign_coords(storm = np.array(l_storms_area))
    xds_TCs_sel.dropna(dim='storm')

    # store storms parameters 
    xds_TCs_sel_params = xr.Dataset(
        {
            'pressure_min':(('storm'), np.array(l_prs_min_in)),
            'pressure_mean':(('storm'), np.array(l_prs_mean_in)),
            'velocity_mean':(('storm'), np.array(l_vel_mean_in)),
            'winds_mean':(('storm'), np.array(l_wind_mean_in)),
            'mean_radius':(('storm'), np.array(l_rad_mean_in)),
            'gamma':(('storm'), np.array(l_gamma)),
            'delta':(('storm'), np.array(l_delta)),
            'category':(('storm'), np.array(l_categ_in)),
            'dmin_date':(('storm'), np.array(l_date_in)),
            #'last_date':(('storm'), np.array(l_date_last)),
            'index_in':(('storm'), np.array(l_ix_in)),
            'index_out':(('storm'), np.array(l_ix_out)),
        },
        coords = {
            'storm':(('storm'), np.arange(cont))
        },
        attrs = {
            'point_lon' : p_lon,
            'point_lat' : p_lat,
            'point_r' : r,
        }
    )

    return xds_TCs_sel, xds_TCs_sel_params


def Extract_Circle(xds_TCs, p_lon, p_lat, r, d_vns):
    '''
    Extracts TCs inside circle - used with NWO or Nakajo databases

    xds_TCs: tropical cyclones track database
        lon, lat, pressure variables
        storm dimension

    circle defined by:
        p_lon, p_lat  -  circle center
        r             -  circle radius (degree)

    d_vns: dictionary to set longitude, latitude, time and pressure varnames

    returns:
        xds_area: selection of xds_TCs inside circle
        xds_inside: contains TCs custom data inside circle
    '''

    # point longitude and latitude
    lonlat_p = np.array([[p_lon, p_lat]])

    # get names of vars: longitude, latitude, pressure and time
    nm_lon = d_vns['longitude']
    nm_lat = d_vns['latitude']
    nm_prs = d_vns['pressure']
    nm_tim = d_vns['time']

    # storms longitude, latitude, pressure and time (if available)
    lon = xds_TCs[nm_lon].values[:]
    lat = xds_TCs[nm_lat].values[:]
    prs = xds_TCs[nm_prs].values[:]
    time = xds_TCs[nm_tim].values[:]
    # get storms inside circle area
    n_storms = xds_TCs.storm.shape[0]
    id_storms = xds_TCs.storm.values
    print(n_storms)
    l_storms_area = []

    # inside parameters holders
    l_prs_min_in = []  # circle minimun pressure
    l_prs_mean_in = []  # circle mean pressure
    l_vel_mean_in = []  # circle mean translation speed
    l_categ_in = []  # circle storm category
    l_date_in = []  # circle date (day)
    l_date_last = []  # last cyclone date
    l_gamma = []  # azimuth
    l_delta = []  # delta

    l_ix_in = []  # historical enters the circle index
    l_ix_out = []  # historical leaves the circle index

    for i_storm in range(n_storms):
        #    for i_storm, st_id in enumerate(id_storms):
        # fix longitude <0 data and skip "one point" tracks
        lon_storm = lon[i_storm]

        if not isinstance(lon_storm, np.ndarray): continue
        lon_storm[lon_storm < 0] = lon_storm[lon_storm < 0] + 360

        # stack storm longitude, latitude
        lonlat_s = np.column_stack(
            (lon_storm, lat[i_storm])
        )

        # index for removing nans
        ix_nonan = ~np.isnan(lonlat_s).any(axis=1)
        lonlat_s = lonlat_s[ix_nonan]

        # calculate geodesic distance (degree)
        geo_dist = []
        for lon_ps, lat_ps in lonlat_s:
            geo_dist.append(GeoDistance(lat_ps, lon_ps, p_lat, p_lon))
        geo_dist = np.asarray(geo_dist)
        #        print(lonlat_s, p_lat, p_lon)
        #        print(geo_dist)

        # find storm inside circle and calculate parameters
        if ((geo_dist < r).any()) & (geo_dist.size > 1):
            #            print(geo_dist.size)
            # storm inside circle
            ix_in = np.where(geo_dist < r)[0][:]

            # storm translation velocity
            geo_dist_ss = []
            for i_row in range(lonlat_s.shape[0] - 1):
                i0_lat, i0_lon = lonlat_s[i_row][1], lonlat_s[i_row][0]
                i1_lat, i1_lon = lonlat_s[i_row + 1][1], lonlat_s[i_row + 1][0]
                geo_dist_ss.append(GeoDistance(i0_lat, i0_lon, i1_lat, i1_lon))
            geo_dist_ss = np.asarray(geo_dist_ss)

            # get delta time in hours (irregular data time delta)
            if isinstance(time[i_storm][0], np.datetime64):
                # round to days
                time[i_storm] = np.array(
                    [np.datetime64(xt, 'h') for xt in time[i_storm]]
                )

                delta_h = np.diff(
                    time[i_storm][~np.isnat(time[i_storm])]
                ).astype('timedelta64[h]').astype(float)

            else:
                # nakajo: time already in hours
                delta_h = np.diff(
                    time[i_storm][~np.isnan(time[i_storm])]
                ).astype(float)

            vel = geo_dist_ss * 111.0 / delta_h  # km/h

            # promediate vel
            velpm = (vel[:-1] + vel[1:]) / 2

            #            print(i_storm, vel.shape, vel[0])
            velpm = np.append(vel[0], velpm)
            velpm = np.append(velpm, vel[-1])

            # calculate azimuth
            lat_in_end, lon_in_end = lonlat_s[ix_in[-1]][1], lonlat_s[ix_in[-1]][0]
            lat_in_ini, lon_in_ini = lonlat_s[ix_in[0]][1], lonlat_s[ix_in[0]][0]
            gamma = GeoAzimuth(lat_in_end, lon_in_end, lat_in_ini, lon_in_ini)
            if gamma < 0.0: gamma += 360

            # calculate delta
            nd = 1000
            st = 2 * np.pi / nd
            ang = np.arange(0, 2 * np.pi + st, st)
            xps = r * np.cos(ang) + p_lat
            yps = r * np.sin(ang) + p_lon
            angle_radius = []
            for x, y in zip(xps, yps):
                angle_radius.append(GeoAzimuth(lat_in_end, lon_in_end, x, y))
            angle_radius = np.asarray(angle_radius)

            im = np.argmin(np.absolute(angle_radius - gamma))
            delta = GeoAzimuth(p_lat, p_lon, xps[im], yps[im])  # (-180, +180)
            if delta < 0.0: delta += 360

            # more parameters
            prs_s_in = prs[i_storm][ix_in]  # pressure
            prs_s_min = np.min(prs_s_in)  # pressure minimun
            prs_s_mean = np.mean(prs_s_in)

            vel_s_in = velpm[ix_in]  # velocity
            vel_s_mean = np.mean(vel_s_in)  # velocity mean

            categ = GetStormCategory(prs_s_min)  # category

            dist_in = geo_dist[ix_in]
            p_dm = np.where((dist_in == np.min(dist_in)))[0]  # closest to point

            time_s_in = time[i_storm][ix_in]  # time
            time_closest = time_s_in[p_dm][0]  # time closest to point

            # filter storms
            # TODO: storms with only one track point inside radius. solve?
            if np.isnan(np.array(prs_s_in)).any() or \
                    (np.array(prs_s_in) <= 860).any() or \
                    gamma == 0.0:
                continue

            # store parameters
            l_storms_area.append(i_storm)
            l_prs_min_in.append(np.array(prs_s_min))
            l_prs_mean_in.append(np.array(prs_s_mean))
            l_vel_mean_in.append(np.array(vel_s_mean))
            l_categ_in.append(np.array(categ))
            l_date_in.append(time_closest)
            l_gamma.append(gamma)
            l_delta.append(delta)

            # store historical indexes inside circle
            l_ix_in.append(ix_in[0])
            l_ix_out.append(ix_in[-1])

            # store last cyclone date too
            l_date_last.append(time[i_storm][ix_nonan][-1])

    #    print(l_storms_area)
    #    print(type(l_storms_area))
    # cut storm dataset to selection
    xds_TCs_sel = xds_TCs.isel(storm=l_storms_area)
    xds_TCs_sel = xds_TCs_sel.assign_coords(storm=np.array(l_storms_area + xds_TCs.storm.values[0]))

    # store storms parameters
    xds_TCs_sel_params = xr.Dataset(
        {
            'pressure_min': (('storm'), np.array(l_prs_min_in)),
            'pressure_mean': (('storm'), np.array(l_prs_mean_in)),
            'velocity_mean': (('storm'), np.array(l_vel_mean_in)),
            'gamma': (('storm'), np.array(l_gamma)),
            'delta': (('storm'), np.array(l_delta)),
            'category': (('storm'), np.array(l_categ_in)),
            'dmin_date': (('storm'), np.array(l_date_in)),
            'last_date': (('storm'), np.array(l_date_last)),
            'index_in': (('storm'), np.array(l_ix_in)),
            'index_out': (('storm'), np.array(l_ix_out)),
        },
        coords={
            'storm': (('storm'), np.array(l_storms_area + xds_TCs.storm.values[0]))
        },
        attrs={
            'point_lon': p_lon,
            'point_lat': p_lat,
            'point_r': r,
        }
    )

    return xds_TCs_sel, xds_TCs_sel_params


def Extract_Circle_list(ls_xds_TCs, p_lon, p_lat, r, d_vns):
    '''
    Extracts TCs inside circle, using a data list, merges a single dataset and
    accomodates attributes, variables...
    '''
    ls_tracks, ls_params = [], []
    for i in range(len(ls_xds_TCs)):
        # extract circle
        xds_TCs_sel, xds_TCs_sel_params = Extract_Circle(ls_xds_TCs[i], p_lon, p_lat, r, d_vns)
        ls_tracks.append(xds_TCs_sel)
        ls_params.append(xds_TCs_sel_params)

    # merge datasets
    xds_TCs_sel = xr.merge(ls_tracks)
    xds_TCs_sel_params = xr.merge(ls_params)

    # add lost attributes
    xds_TCs_sel_params.attrs['point_lon'] = p_lon
    xds_TCs_sel_params.attrs['point_lat'] = p_lat
    xds_TCs_sel_params.attrs['point_r'] = r

    # index_in, index_out and category dtype changed at the merge!!
    xds_TCs_sel_params['index_in'] = xds_TCs_sel_params.index_in.astype('int64')
    xds_TCs_sel_params['index_out'] = xds_TCs_sel_params.index_out.astype('int64')
    xds_TCs_sel_params['category'] = xds_TCs_sel_params.category.astype('int64')

    return xds_TCs_sel, xds_TCs_sel_params


def GetStormCategory(pres_min):
    '''
    Returns storm category (int 5-0)
    '''

    pres_lims = [920, 944, 964, 979, 1000]

    if pres_min <= pres_lims[0]:
        return 5
    elif pres_min <= pres_lims[1]:
        return 4
    elif pres_min <= pres_lims[2]:
        return 3
    elif pres_min <= pres_lims[3]:
        return 2
    elif pres_min <= pres_lims[4]:
        return 1
    else:
        return 0

def SortCategoryCount(np_categ, nocat=9):
    '''
    Sort category change - count matrix
    np_categ = [[category1, category2, count], ...]
    '''

    categs = [0,1,2,3,4,5,9]

    np_categ = np_categ.astype(int)
    np_sort = np.empty((len(categs)*(len(categs)-1),3))
    rc=0
    for c1 in categs[:-1]:
        for c2 in categs:
            p_row = np.where((np_categ[:,0]==c1) & (np_categ[:,1]==c2))
            if p_row[0].size:
                np_sort[rc,:]=[c1,c2,np_categ[p_row,2]]
            else:
                np_sort[rc,:]=[c1,c2,0]

            rc+=1

    return np_sort.astype(int)

def GetCategoryChangeProbs(xds_r1, xds_r2):
    'Calculates category change probabilities from r1 to r2'

    # Get storm category inside both circles
    n_storms = len(xds_r1.storm)
    categ_r1r2 = np.empty((n_storms, 2))
    for i in range(len(xds_r1.storm)):

        # category inside R1
        storm_in_r1 = xds_r1.isel(storm=[i])
        storm_id = storm_in_r1.storm.values[0]
        storm_cat_r1 = storm_in_r1.category

        # category inside R2
        if storm_id in xds_r2.storm.values[:]:
            storm_in_r2 = xds_r2.sel(storm=[storm_id])
            storm_cat_r2 = storm_in_r2.category
        else:
            storm_cat_r2 = 9  # no category

        # store categories
        categ_r1r2[i,:] = [storm_cat_r1, storm_cat_r2]

    # count category changes and sort it
    categ_count = GetUniqueRows(categ_r1r2)
    categ_count = SortCategoryCount(categ_count)

    # calculate probability
    m_count = np.reshape(categ_count[:,2], (6,-1)).T
    m_sum = np.sum(m_count,axis=0)

    probs = m_count.astype(float)/m_sum.astype(float)
    probs_cs = np.cumsum(probs, axis=0)

    # TODO: category_change_sum ??

    # store output using xarray
    xds_categ_cp = xr.Dataset(
        {
            'category_change_count': (('category','category'), m_count[:-1,:]),
            #'category_change_sum': (('category'), m_count[-1,:]),
            'category_change_probs': (('category','category'), probs[:-1,:]),
            'category_nochange_probs': (('category'), probs[-1,:]),
            'category_change_cumsum': (('category','category'), probs_cs[:-1,:]),
        },
        coords = {
            'category': [0,1,2,3,4,5]
        }
    )

    return xds_categ_cp

