#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.colors as mcolors
from matplotlib import cm

def colors_awt():

    # 6 AWT colors
    l_colors_dwt = [
        (155/255.0, 0, 0),
        (1, 0, 0),
        (255/255.0, 216/255.0, 181/255.0),
        (164/255.0, 226/255.0, 231/255.0),
        (0/255.0, 190/255.0, 255/255.0),
        (51/255.0, 0/255.0, 207/255.0),
    ]

    return np.array(l_colors_dwt)

def colors_mjo():
    'custom colors for MJO 25 categories'

    l_named_colors = [
        'lightskyblue', 'deepskyblue', 'royalblue', 'mediumblue',
        'darkblue', 'darkblue', 'darkturquoise', 'turquoise',
        'maroon', 'saddlebrown', 'chocolate', 'gold', 'orange',
        'orangered', 'red', 'firebrick', 'Purple', 'darkorchid',
        'mediumorchid', 'magenta', 'mediumslateblue', 'blueviolet',
        'darkslateblue', 'indigo', 'darkgray',
    ]

    # get rgb colors as numpy array
    np_colors_rgb = np.array(
        [mcolors.to_rgb(c) for c in l_named_colors]
    )

    return np_colors_rgb

def colors_dwt(num_clusters):

    # 42 DWT colors
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

    # get first N colors 
    np_colors_base = np.array(l_colors_dwt)
    np_colors_rgb = np_colors_base[:num_clusters]

    return np_colors_rgb

def colors_interp(num_clusters):

    # generate spectral colormap
    scm = cm.get_cmap('Spectral', num_clusters)

    # use normalize
    mnorm = mcolors.Normalize(vmin=0, vmax=num_clusters)

    # interpolate colors from cmap
    l_colors = []
    for i in range(num_clusters):
        l_colors.append(scm(mnorm(i)))

    # return numpy array
    np_colors_rgb = np.array(l_colors)[:,:-1]

    return np_colors_rgb


def GetClusterColors(num_clusters):
    'Choose colors or Interpolate custom colormap to number of clusters'

    if num_clusters == 6:
        np_colors_rgb = colors_awt()  # Annual Weather Types

    if num_clusters == 25:
        np_colors_rgb = colors_mjo()  # MJO Categories

    elif num_clusters in [36, 42]:
        np_colors_rgb = colors_dwt(num_clusters)  # Daily Weather Types

    else:
        np_colors_rgb = colors_interp(num_clusters)  # interpolate

    return np_colors_rgb
