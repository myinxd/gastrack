# Copyright (C) 2016 Zhixian MA <zxma_sjtu@qq.com>
"""
Tools for simulation

gen_mosaic: make projections of the 3D particles to 2D with respect to the
            direction.
get_peaks: detect peaks in the 2D map and output peaklist
cmp_id: count number of same particles between two different snaps.
"""

import numpy as np


def gen_mosaic(cords, data, step=0.005, direction='z'):
    """
    Make projections and mosaicalize the maps

    Parameters
    ----------
    cords: np.ndarray
        The M*3 matrix that holds the coordinates of the particles.
    data: np.ndarray
        Data to be projected, usually are the masses or density.
    step: double
        length of the mosaic cubic
    direction: str
        Projection direction, default as 'z'.

    Output
    ------
    DensityMat: np.ndarray
        The projected and mosaiced map.
    """
    # Init
    cord_x = cords[:, 0]
    cord_y = cords[:, 1]
    cord_z = cords[:, 2]

    # Cut
    wlim = 2 / 2  # limitation of the box size, 2 Mpc
    # Get indices
    idx_x = (cord_x >= -wlim) * (cord_x <= wlim)
    idx_y = (cord_y >= -wlim) * (cord_y <= wlim)
    idx_z = (cord_z >= -wlim) * (cord_z <= wlim)
    idx = idx_x * idx_y * idx_z

    # Judge direction
    if direction == 'z':
        cord_col = cord_x[idx]
        cord_row = cord_y[idx]
    elif direction == 'x':
        cord_col = cord_y[idx]
        cord_row = cord_z[idx]
    elif direction == 'y':
        cord_col = cord_x[idx]
        cord_row = cord_z[idx]

    data = data[idx]
    # Moscaic and projection
    bins = np.arange(-wlim, wlim + step, step)
    scale = len(bins) - 1
    DensityMat = np.zeros((scale, scale))
    for i in range(scale):
        bin_idx = (cord_col >= bins[i]) * (cord_col <= bins[i + 1])
        temp_row = cord_row[bin_idx]
        temp_data = data[bin_idx]
        for j in range(scale):
            bin_idy = (temp_row >= bins[j]) * (temp_row <= bins[j + 1])
            DensityMat[j, i] = np.sum(temp_data[bin_idy]) / step**2

    return DensityMat


def get_peaks(DensityMat, step=0.005, box_size=2):
    """
    Detect peaks in the map.

    Parameter
    ---------
    DensityMat: np.ndarray
        Two dimensional map, the peak in which shall be detected.
    step: double
        Step between two indices so as to map matrix index to real coordinate.


    Output
    ------
    peaklist: list
        List of the peak
    """
    # Find peak
    peak = np.max(DensityMat)
    idx_row, idx_col = np.where(DensityMat == peak)

    idx_row = np.mean(idx_row)
    idx_col = np.mean(idx_col)

    # Translate
    wlim = box_size / 2
    idx_row = idx_row * step - wlim
    idx_col = idx_col * step - wlim

    # list
    peaklist = [peak, idx_col, idx_row]

    return peaklist


def cmp_id(cls1, cls2, idx1, idx2):
    """Compare same particles between two clusters and output numbers

    Parameters
    ----------
    cls1,cls2: Cluster object
    idx1,idx2: Indices of detected particles in the clusters.

    Output
    ------
    The number of same particles.
    """
    partId1 = cls1.gas_id[idx1]
    partId2 = cls2.gas_id[idx2]
    sameId = np.intersect1d(partId1, partId2)

    return len(sameId)
