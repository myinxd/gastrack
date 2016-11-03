# Copyright (C) 2016 Zhixian MA <zxma_sjtu@qq.com>
# MIT license

"""
Calculate the gas densities with respect to the coordinates provided.

Parameters
----------
posmat: (np.ndarray,unit)
    The position or coordinates of the interested points.
fname: string
    Path of the hdf5 files
field: tuple
    Interestd fields.
"""

import numpy as np
import astro.units as au
import yt
import sys


def calcDensity(fname,posmat,field=('gas','density')):
    # Init
    unit_base = {'UnitLength_in_cm': 3.08568E21,
                 'UnitMass_in_g': 1.989E43,
                 'UnitVelocity_in_cm_per_s': 100000
                }
    # Bound box
    bbox_lim = 1E5 #[kpc]
    bbox = [[-bbox_lim,bbox_lim],
            [-bbox_lim,bbox_lim],
            [-bbox_lim,bbox_lim]]
    # load file and get ds
    try:
        ds = yt.load(fname,unit_base=unit_base,bounding_box=bbox)
    except IOError:
        print('File does not exist.')
        return
    # get field
    ad = ds.all_data()
    cords = np.array(ad['PartType1','Coordinates'])
    if field not in ds.field_list:
        if field in ds.derived_field_list:
            den = ad[field]
        else:
            print("Field does not exist.")
    else:
        den = ad[field]
    # Calc density
    if posmat[1] == 'Mpc':
        pos = posmat[0] * au.Mpc
        pos = pos.to(au.cm)
    elif posmat[1] == 'cm':
        pos = posmat[0] * au.cm
    numpoints = pos.shape[0]
    posdens = np.zeros((numpoints,))
    for i in range(numpoints):
        dist = (cords[:,0]-pos[i,0])**2 + (cords[:,1]-pos[i,1])**2
        dist = np.sqrt(dist)
        idx = np.where(dist == dist.min())[0]
        posdens[i] = den[idx]

    return posdens

def main(argv):
    """
    A main function to run this script at commind lines
    """
    fname = argv[0]
    pos_in_str = argv[1:-1].split()
    unit = argv[-1:]
    num = len(pos_in_str)//2
    # Get point matrix
    posmat = np.zeros((num,2))
    for i in range(num):
        posmat[i,:] = pos_in_str[i*2:i*2+1]
    # Calc density
    posdens = calcDensity(fname,(posmat,unit))
    # output
    for i in range(num):
        print("coord: (%f, %f) Mpc, Density: %f g/cm3" %
              (posmat[i,:],posdens[i]))


if __name__== "__main__":
    main(sys.argv)
