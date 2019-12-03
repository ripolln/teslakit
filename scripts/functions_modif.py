
import matplotlib.pyplot as plt
from teslakit.plotting.config import _faspect, _fsize, _fdpi
import numpy as np

def Plot_DWTs_Mean_Anom_modif(xds_KMA, xds_var, kind='mean', p_export=None):
    '''
    Plot Daily Weather Types (bmus mean)
    kind - mean/anom
    '''

    bmus = xds_KMA['sorted_bmus'].values[:]
    n_clusters = len(xds_KMA.n_clusters.values[:])

    # plot figure
    fig, ax = plt.subplots()

    values = []
    for ic in range(n_clusters):

        if kind=='mean':
            # data mean
            it = np.where(bmus==ic)[0][:]
            c_plot = xds_var.isel(time=it).mean(dim='time')

        elif kind=='anom':
            # data anomally
            it = np.where(bmus==ic)[0][:]
            t_mean = xds_var.mean(dim='time')
            c_mean = xds_var.isel(time=it).mean(dim='time')
            c_plot = c_mean - t_mean

        values = np.append(values, c_plot.values)

    values = np.reshape(values,(7,6))
    values = np.flipud(values)

    plt.pcolor(values)
    plt.colorbar()
    plt.clim(0,8)
    ax.tick_params(axis='both', which='both', length=0)
    ax.set_xticklabels('')
    ax.set_yticklabels('')

    plt.title('mean Area below TWL max (m*h)')

    # show / export
    if not p_export:
        plt.show()
    else:
        fig.savefig(p_export, dpi=_fdpi)
        plt.close()


def Plot_Probs_WT_MJO(series_1, series_2, n_clusters_1, n_clusters_2, ttl='',
                     wt_colors=False, show=True, rows=5, cols=5):
    '''
    Plot WTs_1 / WTs_2 probabilities

    both categories series should start at 0
    '''

    # set of daily weather types
    set_2 = np.arange(n_clusters_2)

    # dailt weather types matrix rows and cols
    n_rows, n_cols = GetBestRowsCols(n_clusters_2)

    # get cluster colors
    cs_wt =  GetClusterColors(n_clusters_1)

    # plot figure
    fig = plt.figure(figsize=(_faspect*_fsize, _fsize))
    gs = gridspec.GridSpec(rows, cols, wspace=0.10, hspace=0.15)

    xpos=0
    ypos=-1
    for ic in range(n_clusters_1):

        # select DWT bmus at current AWT indexes
        index_1 = np.where(series_1==ic)[0][:]
        sel_2 = series_2[index_1]

        # get DWT cluster probabilities
        cps = ClusterProbabilities(sel_2, set_2)
        C_T = np.reshape(cps, (n_rows, n_cols))
        print(ic+1)
        print(C_T)
        print()
        print()

        # axis colors
        if wt_colors:
            caxis = cs_wt[ic]
        else:
            caxis = 'black'

        # contador
        ypos +=1
        if ypos >= cols:
            xpos +=1
            ypos=0

        # plot axes
        ax = plt.subplot(gs[xpos, ypos])
        axplot_WT_Probs(
            ax, C_T,
            ttl = 'WT {0}'.format(ic+1),
            cmap = 'Reds', caxis = caxis,
        )
        ax.set_aspect('equal')

    # add fig title
    fig.suptitle(ttl, fontsize=14, fontweight='bold')
    plt.tight_layout()

    # show and return figure
    if show: plt.show()
    return fig
