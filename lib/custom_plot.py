#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.interpolate import interp1d
import calendar
from datetime import datetime, timedelta

from lib.util.terminal import printProgressBar as pb

def Plot_PredictorEOFs(xds_PCA, n_plot):
    '''
    Plot EOFs
    '''
    # TODO: DOC

    # PCA data
    variance = xds_PCA['variance'].values
    EOFs = np.transpose(xds_PCA['EOFs'].values)
    PCs = np.transpose(xds_PCA['PCs'].values)

    years = xds_PCA['_years'].values
    lon = xds_PCA['_longitude'].values
    len_x = len(lon)

    m1 = xds_PCA.attrs['m1']
    m2 = xds_PCA.attrs['m2']
    l_months = [calendar.month_name[x] for x in range(1,13)]
    ylbl = l_months[m1-1:] + l_months[:m2]

    # percentage of variance each field explains
    n_percent = variance / np.sum(variance)

    for it in range(n_plot):

        # map of the spatial field
        spatial_fields = EOFs[:,it]*np.sqrt(variance[it])

        # reshape from vector to matrix with separated months
        C = np.reshape(spatial_fields[:len_x*12], (12, len_x)).transpose()

        # eof cmap
        ax1 = plt.subplot2grid((6, 6), (0, 0), colspan=6, rowspan=4)
        plt.pcolormesh(np.transpose(C), cmap='RdBu', shading='gouraud')
        plt.clim(-1,1)
        plt.title('EOF #{0}  ---  {1:.2f}%'.format(it+1,n_percent[it]*100))
        ax1.set_xticklabels([str(x) for x in lon])
        ax1.set_yticklabels(ylbl)

        # time series
        ax2 = plt.subplot2grid((6, 6), (5, 0), colspan=6, rowspan=2)
        plt.plot(years, PCs[it,:]/np.sqrt(variance[it]))
        plt.xlim(years[0], years[-1])

        # SHOW
        plt.show()

def Plot_MJOphases(rmm1, rmm2, phase):
    'Plot MJO data separated by phase'

    # parameters for custom plot
    size_points = 0.2
    size_lines = 0.8
    l_colors_phase = [
        (1, 0, 0),
        (0.6602, 0.6602, 0.6602),
        (1.0, 0.4961, 0.3125),
        (0, 1, 0),
        (0.2539, 0.4102, 0.8789),
        (0, 1, 1),
        (1, 0.8398, 0),
        (0.2930, 0, 0.5078)]
    color_lines_1 = (0.4102, 0.4102, 0.4102)


    # plot data
    plt.figure(1)
    ax = plt.subplot(111)
    ax.scatter(rmm1, rmm2, c='b', s=size_points)

    # plot data by phases
    for i in range(1,9):
        ax.scatter(
            rmm1.where(phase==i),
            rmm2.where(phase==i),
            c=l_colors_phase[i-1],
            s=size_points)

    # plot sectors
    ax.plot([-4,4],[-4,4], color='k', linewidth=size_lines)
    ax.plot([-4,4],[4,-4], color='k', linewidth=size_lines)
    ax.plot([-4,4],[0,0],  color='k', linewidth=size_lines)
    ax.plot([0,0], [-4,4], color='k', linewidth=size_lines)

    # axis
    plt.xlim(-4, 4)
    plt.ylim(-4, 4)
    plt.xlabel('RMM1')
    plt.ylabel('RMM2')
    ax.set_aspect('equal')

    # show
    plt.show()

def Plot_MJOCategories(rmm1, rmm2, categ):
    'Plot MJO data separated by 25 categories'

    # parameters for custom plot
    size_lines = 0.8
    color_lines_1 = (0.4102, 0.4102, 0.4102)
    # TODO: COLORES PARA 25 CATEGORIAS, NO PARA N
    l_colors_categ = [
        (0.527343750000000, 0.804687500000000, 0.979166666666667),
        (0, 0.746093750000000, 1),
        (0.253906250000000, 0.410156250000000, 0.878906250000000),
        (0, 0, 0.800781250000000),
        (0, 0, 0.542968750000000),
        (0.273437500000000, 0.507812500000000, 0.703125000000000),
        (0, 0.804687500000000, 0.816406250000000),
        (0.250000000000000, 0.875000000000000, 0.812500000000000),
        (0.500000000000000, 0, 0),
        (0.542968750000000, 0.269531250000000, 0.0742187500000000),
        (0.820312500000000, 0.410156250000000, 0.117187500000000),
        (1, 0.839843750000000, 0),
        (1, 0.644531250000000, 0),
        (1, 0.269531250000000, 0),
        (1, 0, 0),
        (0.695312500000000, 0.132812500000000, 0.132812500000000),
        (0.500000000000000, 0, 0.500000000000000),
        (0.597656250000000, 0.195312500000000, 0.796875000000000),
        (0.726562500000000, 0.332031250000000, 0.824218750000000),
        (1, 0, 1),
        (0.480468750000000, 0.406250000000000, 0.929687500000000),
        (0.539062500000000, 0.167968750000000, 0.882812500000000),
        (0.281250000000000, 0.238281250000000, 0.542968750000000),
        (0.292968750000000, 0, 0.507812500000000),
        (0.660156250000000, 0.660156250000000, 0.660156250000000),
    ]

    # plot figure
    plt.figure(1)
    ax = plt.subplot(111)

    # plot sectors
    ax.plot([-4,4],[-4,4], color='k', linewidth=size_lines, zorder=9)
    ax.plot([-4,4],[4,-4], color='k', linewidth=size_lines, zorder=9)
    ax.plot([-4,4],[0,0],  color='k', linewidth=size_lines, zorder=9)
    ax.plot([0,0], [-4,4], color='k', linewidth=size_lines, zorder=9)

    # plot circles
    R = [1, 1.5, 2.5]

    for rr in R:
        ax.add_patch(
            patches.Circle(
                (0,0),
                rr,
                color='k',
                linewidth=size_lines,
                fill=False,
                zorder=9)
        )
    ax.add_patch(
        patches.Circle((0,0),R[0],fc='w',fill=True, zorder=10))

    # plot data by categories
    for i in range(1,25):
        if i>8: size_points = 0.2
        else: size_points = 1.7

        ax.scatter(
            rmm1.where(categ==i),
            rmm2.where(categ==i),
            c=l_colors_categ[i-1],
            s=size_points)
    ax.scatter(
        rmm1.where(categ==25),
        rmm2.where(categ==25),
        c=l_colors_categ[i-1],
        s=0.2,
    zorder=11)

    # axis
    plt.xlim(-4, 4)
    plt.ylim(-4, 4)
    plt.xlabel('RMM1')
    plt.ylabel('RMM2')
    ax.set_aspect('equal')

    # show
    plt.show()

def Plot_ARL_PerpYear(bmus_values, bmus_dates, num_clusters, num_sims):
    'Plots ARL simulation output in a perpetual_year stacked bar chart'

    # parameters for custom plot
    l_colors_dwt = [
        (1.0000, 0.1344, 0.0021),
        (1.0000, 0.2669, 0.0022),
        (1.0000, 0.5317, 0.0024),
        (1.0000, 0.6641, 0.0025),
        (1.0000, 0.9287, 0.0028),
        (0.9430, 1.0000, 0.0029),
        (0.6785, 1.0000, 0.0031),
        (0.5463, 1.0000, 0.0032),
        (0.2821, 1.0000, 0.0035),
        (0.1500, 1.0000, 0.0036),
        (0.0038, 1.0000, 0.1217),
        (0.0039, 1.0000, 0.2539),
        (0.0039, 1.0000, 0.4901),
        (0.0039, 1.0000, 0.6082),
        (0.0039, 1.0000, 0.8444),
        (0.0039, 1.0000, 0.9625),
        (0.0039, 0.8052, 1.0000),
        (0.0039, 0.6872, 1.0000),
        (0.0040, 0.4510, 1.0000),
        (0.0040, 0.3329, 1.0000),
        (0.0040, 0.0967, 1.0000),
        (0.1474, 0.0040, 1.0000),
        (0.2655, 0.0040, 1.0000),
        (0.5017, 0.0040, 1.0000),
        (0.6198, 0.0040, 1.0000),
        (0.7965, 0.0040, 1.0000),
        (0.8848, 0.0040, 1.0000),
        (1.0000, 0.0040, 0.9424),
        (1.0000, 0.0040, 0.8541),
        (1.0000, 0.0040, 0.6774),
        (1.0000, 0.0040, 0.5890),
        (1.0000, 0.0040, 0.4124),
        (1.0000, 0.0040, 0.3240),
        (1.0000, 0.0040, 0.1473),
        (0.9190, 0.1564, 0.2476),
        (0.7529, 0.3782, 0.4051),
        (0.6699, 0.4477, 0.4584),
        (0.5200, 0.5200, 0.5200),
        (0.4595, 0.4595, 0.4595),
        (0.4100, 0.4100, 0.4100),
        (0.3706, 0.3706, 0.3706),
        (0.2000, 0.2000, 0.2000),
        (     0, 0, 0),
    ]

    # TODO: GUARDAR EL INTERPOLADOR EN UNA FUNCION
    # interpolate colors to num cluster
    np_colors_base = np.array(l_colors_dwt)
    x = np.arange(np_colors_base.shape[0])
    itp = interp1d(x, np_colors_base, axis=0, kind='linear')

    xi = np.arange(num_clusters)
    np_colors_int =  itp(xi)



    # generate perpetual year list
    dp1 = datetime(1981,1,1)
    dp2 = datetime(1981,12,31)
    list_pyear = [dp1 + timedelta(days=i) for i in range((dp2-dp1).days+1)]

    # generate aux arrays
    m_plot = np.zeros((num_clusters, len(list_pyear))) * np.nan
    bmus_dates_months = np.array([d.month for d in bmus_dates])
    bmus_dates_days = np.array([d.day for d in bmus_dates])

    # sort data
    for i, dpy in enumerate(list_pyear):
        _, s = np.where(
            [(bmus_dates_months == dpy.month) & (bmus_dates_days == dpy.day)]
        )
        b = bmus_values[s]

        for j in range(num_clusters):
            _, bb = np.where([(j+1 == b)])

            m_plot[j,i] = float(len(bb))/len(s)


    # plot figure
    plt.figure(1)
    ax = plt.subplot(111)

    bottom_val = np.zeros(m_plot[1,:].shape)
    for r in range(num_clusters):
        row_val = m_plot[r,:]
        plt.bar(
            range(365), row_val, bottom=bottom_val,
            width=1, color = np_colors_int[r]
               )

        # store bottom
        bottom_val += row_val

        # progressbar
        pb(r, num_clusters,
            prefix = 'Plotting',
            suffix = 'Complete', length = 50)

    # axis
    plt.xlim(1, 365)
    plt.ylim(0, 1)
    plt.xlabel('Perpetual year')
    plt.ylabel('')

    # show
    plt.show()



