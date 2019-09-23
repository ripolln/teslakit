
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import numpy as np
import copy
from datetime import datetime

# pip
from cftime._cftime import DatetimeGregorian
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
from matplotlib import cm
import pandas as pd

# teslakit
from .custom_colors import colors_dwt
from ..util.operations import GetBestRowsCols
from ..custom_dateutils import npdt64todatetime as n2d
from ..kma import ClusterProbabilities
from .awt import axplot_PCs_3D_allWTs

# import constants
from .config import _faspect, _fsize, _fdpi

def axplot_EOF(ax, EOF_value, lon, lat, ttl=''):
    'axes plot EOFs 2d map'

    cmap = cm.get_cmap('RdBu_r')

    # EOF pcolormesh 
    ax.pcolormesh(
        lon, lat, np.transpose(EOF_value),
        cmap=cmap, shading='gouraud',
        clim=(-1,1),
    )

    # axis and title
    ax.set_title(
        ttl,
        {'fontsize': 10, 'fontweight':'bold'}
    )
    ax.tick_params(axis='both', which='major', labelsize=8)

def axplot_EOF_evolution(ax, time, EOF_evol):
    'axes plot EOFs evolution'

    # date axis locator
    yloc1 = mdates.YearLocator(1)
    yfmt = mdates.DateFormatter('%Y')

    # convert to datetime
    dtime = [n2d(t) for t in time]

    # plot EOF evolution 
    ax.plot(
        dtime, EOF_evol,
        linestyle='-', linewidth=0.5, color='black',
    )

    # configure axis
    ax.set_xlim(time[0], time[-1])
    ax.xaxis.set_major_locator(yloc1)
    ax.xaxis.set_major_formatter(yfmt)
    ax.grid(True, which='both', axis='x', linestyle='--', color='grey')
    ax.tick_params(axis='both', which='major', labelsize=8)

def axplot_DWT(ax, dwt, vmin, vmax, wt_color):
    'axes plot EOFs 2d map'

    cmap = copy.deepcopy(cm.get_cmap('RdBu_r'))
    cmap.set_bad(color = wt_color, alpha=0.2)

    # EOF pcolormesh 
    pc = ax.pcolormesh(
        dwt,
        cmap = cmap, shading = 'gouraud',
        clim = (vmin, vmax),
    )

    # customize axis
    ax.set_xticks([])
    ax.set_yticks([])

    return pc

def axplot_DWT_Probs(ax, dwt_probs,
                     ttl = '', vmin = 0, vmax = 0.1,
                     cmap = 'Blues'):
    'axes plot DWT cluster probabilities'

    # clsuter transition plot
    pc = ax.pcolor(
        np.flipud(dwt_probs),
        cmap=cmap, vmin=vmin, vmax=vmax,
        edgecolors='k',
    )

    # customize axes
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(ttl, {'fontsize':10, 'fontweight':'bold'})

    return pc

def axplot_DWT_Hist(ax, bmus, n_clusters, ttl=''):
    'axes plot DWT cluster count histogram'

    # clsuter transition plot
    ax.hist(
        bmus,
        bins = np.arange(1, n_clusters+2),
        edgecolor='k'
    )

    # customize axes
    #ax.grid('y')

    ax.set_xticks(np.arange(1,n_clusters+1)+0.5)
    ax.set_xticklabels(np.arange(1,n_clusters+1))
    ax.set_xlim([1, n_clusters+1])
    ax.tick_params(axis='both', which='major', labelsize=6)

    ax.set_title(ttl, {'fontsize':10, 'fontweight':'bold'})


def Plot_EOFs_EstelaPred(xds_PCA, n_plot, p_export=None):
    '''
    Plot annual EOFs for 3D predictors

    xds_PCA:
        (n_components, n_components) PCs
        (n_components, n_features) EOFs
        (n_components, ) variance

        (n_lon, ) pred_lon: predictor longitude values
        (n_lat, ) pred_lat: predictor latitude values
        (n_time, ) pred_time: predictor time values

        method: gradient + estela

    n_plot: number of EOFs plotted

    show plot or saves figure to p_export
    '''

    # PCA data
    variance = xds_PCA['variance'].values[:]
    EOFs = np.transpose(xds_PCA['EOFs'].values[:])
    PCs = np.transpose(xds_PCA['PCs'].values[:])
    data_pos = xds_PCA['pred_data_pos'].values[:]  # for handling nans
    pca_time = xds_PCA['pred_time'].values[:]
    pred_name = xds_PCA.attrs['pred_name']

    # PCA lat lon metadata
    lon = xds_PCA['pred_lon'].values
    lat = xds_PCA['pred_lat'].values

    # percentage of variance each field explains
    n_percent = variance / np.sum(variance)

    for it in range(n_plot):

        # get vargrd 
        var_grd_1d = EOFs[:,it] * np.sqrt(variance[it])

        # insert nans in data
        base = np.nan * np.ones(data_pos.shape)
        base[data_pos] = var_grd_1d


        var = base[:int(len(base)/2)]
        grd = base[int(len(base)/2):]

        # reshape data to grid
        C1 = np.reshape(var, (len(lon), len(lat)))
        C2 = np.reshape(grd, (len(lon), len(lat)))

        # figure
        fig = plt.figure(figsize=(_faspect*_fsize, _fsize))

        # layout
        gs = gridspec.GridSpec(4, 4, wspace=0.10, hspace=0.2)

        ax_EOF_1 = plt.subplot(gs[:3, :2])
        ax_EOF_2 = plt.subplot(gs[:3, 2:])
        ax_evol = plt.subplot(gs[3, :])

        # EOF pcolormesh (SLP and GRADIENT)
        axplot_EOF(ax_EOF_1, C1, lon, lat, ttl = pred_name)
        axplot_EOF(ax_EOF_2, C2, lon, lat, ttl = 'GRADIENT')

        # time series EOF evolution
        evol =  PCs[it,:]/np.sqrt(variance[it])
        axplot_EOF_evolution(ax_evol, pca_time, evol)

        # figure title
        ttl = 'EOF #{0}  ---  {1:.2f}%'.format(it+1, n_percent[it]*100)
        fig.suptitle(ttl, fontsize=14, fontweight='bold')

        # show / export
        if not p_export:
            plt.show()

        else:
            if not op.isdir(p_export):
                os.makedirs(p_export)
            p_expi = op.join(p_export, 'EOFs_{0}_{1}.png'.format(pred_name, it+1))
            fig.savefig(p_expi, dpi=_fdpi)
            plt.close()

def Plot_ESTELA(pnt_lon, pnt_lat, estela_D, p_export=None):
    'Plots ESTELA days at world map '

    try:
        from mpl_toolkits.basemap import Basemap
    except:
        print('basemap module required.')
        return

    # estela data
    estela_lon = estela_D.longitude.values[:]
    estela_lat = estela_D.latitude.values[:]
    estela_val = estela_D.values[:]
    estela_max = np.ceil(estela_D.max().values)

    # figure
    fig, ax = plt.subplots(1, figsize=(_faspect*_fsize, _fsize))

    # setup mercator map projection.
    m = Basemap(
        #llcrnrlon = lon1, llcrnrlat = lat1,
        #urcrnrlon = lon2, urcrnrlat = lat2,
        resolution = 'l', projection = 'cyl',
        lat_0 = pnt_lat, lon_0 = pnt_lon,
        area_thresh = 0.01,
    )
    m.drawcoastlines()
    m.fillcontinents(color = 'silver')
    m.drawmapboundary(fill_color = 'lightcyan')
    m.drawparallels(np.arange(-90, 90, 20), labels = [1,1,0,0])
    m.drawmeridians(np.arange(0, 360, 20), labels = [0,0,0,1])

    # plot estela
    pc = ax.pcolormesh(
        estela_lon, estela_lat, estela_val,
        cmap='jet_r', shading='gouraud',
        clim=(0, estela_max),
    )
    cb = m.colorbar(pc, location='bottom')
    cb.set_ticks(range(int(estela_max)))
    cb.set_label('days')

    # plot point
    ax.plot(pnt_lon, pnt_lat, 'ok')

    # show / export
    if not p_export:
        plt.show()
    else:
        fig.savefig(p_export, dpi=_fdpi)
        plt.close()

def Plot_DWTs_Mean(xds_KMA, xds_var, p_export=None):
    '''
    Plot Daily Weather Types (bmus mean)
    '''

    bmus = xds_KMA['sorted_bmus'].values[:]
    n_clusters = len(xds_KMA.n_clusters.values[:])

    var_max = np.max(xds_var.values)
    var_min = np.min(xds_var.values)
    scale = 1/100.0  # scale from Pa to mbar

    # get number of rows and cols for gridplot 
    n_rows, n_cols = GetBestRowsCols(n_clusters)

    # get cluster colors
    cs_dwt = colors_dwt(n_clusters)

    # plot figure
    fig = plt.figure(figsize=(_faspect*_fsize, _fsize))

    gs = gridspec.GridSpec(n_rows, n_cols, wspace=0.0, hspace=0.0)
    gr = 0
    gc = 0

    for ic in range(n_clusters):
        # data mean
        it = np.where(bmus==ic)[0][:]
        c_mean = xds_var.isel(time=it).mean(dim='time')

        # convert input units
        c_mean_s = np.multiply(c_mean, scale)

        # dwt color
        clr = cs_dwt[ic]

        # axes plot
        ax = plt.subplot(gs[gr, gc])
        pc = axplot_DWT(
            ax, np.flipud(c_mean_s),
            vmin = var_min, vmax = var_max, wt_color = clr
        )

        # get lower positions
        if gr==n_rows-1 and gc==0:
            pax_l = ax.get_position()
        elif gr==n_rows-1 and gc==n_cols-1:
            pax_r = ax.get_position()

        # counter
        gc += 1
        if gc >= n_cols:
            gc = 0
            gr += 1

    # add a colorbar        
    cbar_ax = fig.add_axes([pax_l.x0, pax_l.y0-0.05, pax_r.x1 - pax_l.x0, 0.02])
    cb = fig.colorbar(pc, cax=cbar_ax, orientation='horizontal')
    cb.set_label('Pressure (mbar)')

    # show / export
    if not p_export:
        plt.show()
    else:
        fig.savefig(p_export, dpi=_fdpi)
        plt.close()

def ClusterProbs_Month(bmus, time, wt_set, month_ix):
    'Returns Cluster probs by month_ix'

    # TODO: update custom_dateutils library
    # get months
    if isinstance(time[0], DatetimeGregorian):
        tpd_month = np.asarray([t.month for t in time])

    else:
        tpd = pd.DatetimeIndex(time)
        tpd_month = tpd.month

    if isinstance(month_ix, list):

        # get each month indexes
        l_ix = []
        for m_ix in month_ix:
            ixs = np.where(tpd_month==m_ix)[0]
            l_ix.append(ixs)

        # get all indexes     
        ix = np.unique(np.concatenate(tuple(l_ix)))

    else:
        ix = np.where(tpd_month==month_ix)[0]

    bmus_sel = bmus[ix]

    return ClusterProbabilities(bmus_sel, wt_set)

def Plot_DWTs_Probs(bmus, bmus_time, n_clusters, p_export=None):
    '''
    Plot Daily Weather Types bmus probabilities
    '''

    wt_set = np.arange(n_clusters) + 1

    # best rows cols combination
    n_rows, n_cols = GetBestRowsCols(n_clusters)

    # figure
    fig = plt.figure(figsize=(_faspect*_fsize, _fsize))

    # layout
    gs = gridspec.GridSpec(4, 7, wspace=0.10, hspace=0.25)

    # list all plots params
    l_months = [
        (1, 'January',   gs[1,3]),
        (2, 'February',  gs[2,3]),
        (3, 'March',     gs[0,4]),
        (4, 'April',     gs[1,4]),
        (5, 'May',       gs[2,4]),
        (6, 'June',      gs[0,5]),
        (7, 'July',      gs[1,5]),
        (8, 'August',    gs[2,5]),
        (9, 'September', gs[0,6]),
        (10, 'October',  gs[1,6]),
        (11, 'November', gs[2,6]),
        (12, 'December', gs[0,3]),
    ]

    l_3months = [
        ([12, 1, 2],  'DJF', gs[3,3]),
        ([3, 4, 5],   'MAM', gs[3,4]),
        ([6, 7, 8],   'JJA', gs[3,5]),
        ([9, 10, 11], 'SON', gs[3,6]),
    ]

    # plot total probabilities
    c_T = ClusterProbabilities(bmus, wt_set)
    C_T = np.reshape(c_T, (n_rows, n_cols))

    ax_probs_T = plt.subplot(gs[:2, :2])
    pc = axplot_DWT_Probs(ax_probs_T, C_T, ttl = 'DWT Probabilities')

    # plot counts histogram
    ax_hist = plt.subplot(gs[2:, :3])
    axplot_DWT_Hist(ax_hist, bmus, n_clusters, ttl = 'DWT Counts')

    # plot probabilities by month
    for m_ix, m_name, m_gs in l_months:

        # get probs matrix
        c_M = ClusterProbs_Month(bmus, bmus_time, wt_set, m_ix)
        C_M = np.reshape(c_M, (n_rows, n_cols))

        # plot axes
        ax_M = plt.subplot(m_gs)
        axplot_DWT_Probs(ax_M, C_M, ttl = m_name)

    # plot probabilities by 3 month sets
    for m_ix, m_name, m_gs in l_3months:

        # get probs matrix
        c_M = ClusterProbs_Month(bmus, bmus_time, wt_set, m_ix)
        C_M = np.reshape(c_M, (n_rows, n_cols))

        # plot axes
        ax_M = plt.subplot(m_gs)
        axplot_DWT_Probs(ax_M, C_M, ttl = m_name, cmap='Greens')

    # add custom colorbar
    pp = ax_probs_T.get_position()
    cbar_ax = fig.add_axes([pp.x1+0.02, pp.y0, 0.02, pp.y1 - pp.y0])
    cb = fig.colorbar(pc, cax=cbar_ax, cmap='Blues')
    cb.ax.tick_params(labelsize=8)

    # show / export
    if not p_export:
        plt.show()
    else:
        fig.savefig(p_export, dpi=_fdpi)
        plt.close()

def Plot_DWT_PCs_3D(d_PCs, n_clusters, p_export=None):
    '''
    Plot Annual Weather Types PCs fit - rnd comparison (3D)
    '''

    from mpl_toolkits.mplot3d import Axes3D

    # get cluster colors
    cs_awt = colors_dwt(n_clusters)

    # figure
    fig = plt.figure(figsize=(_faspect*_fsize, _fsize*2/1.66))
    gs = gridspec.GridSpec(1, 1, wspace=0.10, hspace=0.35)
    ax_pcs = plt.subplot(gs[0, 0], projection='3d')

    # Plot PCs (3D)
    axplot_PCs_3D_allWTs(ax_pcs, d_PCs,  cs_awt, ttl='DWT PCs')

    # show / export
    if not p_export:
        plt.show()
    else:
        fig.savefig(p_export, dpi=_fdpi)
        plt.close()


# TODO: not updated functions
def Plot_PCvsPC(xds_PC123, text=[], p_export = None):
    '''
    Plot PC1 vs PC2 vs PC3

    xds_PD123
        (dim,) PC1
        (dim,) PC2
        (dim,) PC3

        (dim,) text

    show plot or saves figure to p_export
    '''

    # get data
    pc1_val = xds_PC123.PC1.values
    pc2_val = xds_PC123.PC2.values
    pc3_val = xds_PC123.PC3.values

    # delta axis
    pc1_d = np.max([np.absolute(np.max(pc1_val)), np.absolute(np.min(pc1_val))])
    pc2_d = np.max([np.absolute(np.max(pc2_val)), np.absolute(np.min(pc2_val))])
    pc3_d = np.max([np.absolute(np.max(pc3_val)), np.absolute(np.min(pc3_val))])

    pf = 1.05
    pc1_d = pc1_d * pf
    pc2_d = pc2_d * pf
    pc3_d = pc3_d * pf

    # create figure
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, figsize=(_faspect*_fsize, _fsize))

    ax1.plot(pc2_val, pc1_val, '.r')
    ax2.plot(pc3_val, pc1_val, '.r')
    ax4.plot(pc3_val, pc2_val, '.r')
    ax3.remove()

    # text
    for p1,p2,p3,t in zip(pc1_val,pc2_val,pc3_val,text):
        ax1.text(p2,p1,t)
        ax2.text(p3,p1,t)
        ax4.text(p3,p2,t)

    # labels and customize
    fw = 'bold'
    ax1.set_xlabel('PC2', fontweight=fw)
    ax1.set_ylabel('PC1', fontweight=fw)
    ax2.set_xlabel('PC3', fontweight=fw)
    ax2.set_ylabel('PC1', fontweight=fw)
    ax4.set_xlabel('PC3', fontweight=fw)
    ax4.set_ylabel('PC2', fontweight=fw)

    ax1.set_xlim(-pc2_d, pc2_d)
    ax1.set_ylim(-pc1_d, pc1_d)
    ax2.set_xlim(-pc3_d, pc3_d)
    ax2.set_ylim(-pc1_d, pc1_d)
    ax4.set_xlim(-pc3_d, pc3_d)
    ax4.set_ylim(-pc2_d, pc2_d)

    lc = 'k'
    lls = '--'
    for ax in [ax1, ax2, ax4]:
        ax.axhline(y=0, color=lc, linestyle=lls)
        ax.axvline(x=0, color=lc, linestyle=lls)

    # show / export
    if not p_export:
        plt.show()

    else:
        fig.savefig(p_export, dpi=96)
        plt.close()

