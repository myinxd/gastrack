#!/usr/bin/python3
# Copyright (C) 2016 Zhixian MA <zxma_sjtu@qq.com>
# Version 2.0
"""
A gas tracker to detect halo and gas peaks of the two merging galaxies.
The `hdf5` files of simulated result generated from Gadget are as the
input, and two images of the **halo** and **gas** are the outputs. Gas
peaks, i.e., the particles with highest local mass densities are marked
on the gas image.

Parameters
----------
input_dir: string
    The folder holding those simulated results.
output_dir: string
    The folder to save ouput results.
width: tuple
    Width of the output image, e,g. (2,'Mpc').

References
----------
[1] User-guide of gadget-2
    http://www.gadgetcode.org/
[2] yt documents
    http://yt-project.org/doc/index.html
[3] Matplotlib Axes3D
    http://matplotlib.org/mpl_toolkits/mplot3d/api.html?
    highlight=scatter3d#mpl_toolkits.mplot3d.axes3d.Axes3D.scatter3D
"""

import os
import re
import sys
import numpy as np
import yt
from yt.visualization.fixed_resolution import FixedResolutionBuffer
from scipy.ndimage import imread
import matplotlib.pyplot as plt 
from mpl_toolkit.mplot3d import Axes3D

import partsplit


class GasTrack:

    """
    The gas tracker class.

    Methods
    -------
    get_halo_map:
        Get halo map by the yt YTQuadTreeProj and FixedResolutionBuffer
        objects.
    get_gas_map:
        Generate gas map with peaks marked.
    locate_peak:
        Locate peaks in the halo map.
    load_hdf5:
        Load the hdf5 file

    References
    ----------
    [1] YTQuadTreeProj
        http://yt-project.org/doc/reference/api/generated/yt.data_objects.
        construction_data_containers.YTQuadTreeProj.html?highlight=
        ytquadtreeproj
    [2] FixedResolutionBuffer
        http://yt-project.org/doc/reference/api/generated/yt.visualization.
        fixed_resolution.FixedResolutionBuffer.html?highlight=
        fixedresolutionbuffer
    """

    def __init__(self, input_dir, output_dir, width=(0.8, 'Mpc'),
                 buffsize=(3200, 3200),cmaps=('ocean','Accent'),
                 smooth_width=3):
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.save = True
        # Projection axis
        self.axis = 'z'
        # center
        self.center = [0.0,0.0,0.0]
        # width
        self.width = width
        # buff_size
        self.buff_size = buffsize
        # image size
        self.imgsize = (800, 800)
        # colormaps
        self.cmaps = cmaps
        # smooth_width
        self.smooth = smooth_width
        # unit_base
        self.unit_base = {'UnitLength_in_cm': 3.08568E21,
                          'UnitMass_in_g': 1.989E43,
                          'UnitVelocity_in_cm_per_s': 100000}
        # width limit
        self.wlim = width[0]/2
        # Bound box
        bbox_lim = 1E5  # [kpc]
        self.bbox = [[-bbox_lim, bbox_lim],
                     [-bbox_lim, bbox_lim],
                     [-bbox_lim, bbox_lim]]
        # fields
        self.halo_field = ('deposit', 'PartType1_density')
        self.gas_field = ('PartType0', 'Density')

    def load_hdf5(self, fname):
        """Load the hdf5 files

        Parameter
        ---------
        fname: str
            Filename.

        Return
        ------
        ds: yt.frontends.gadget.data_structures.GadgetHDF5Dataset
            The yt GadgetHDF5Dataset object contains fields we are
            interested.
        """
        fpath = os.path.join(self.input_dir, fname)
        ds = yt.load(fpath, unit_base=self.unit_base,
                     bounding_box=self.bbox)
        return ds

    def get_window_parameters(self, ds):
        """Get bounds of the axes."""

        width = ds.coordinates.sanitize_width(self.axis, self.width, None)
        center, display_center = ds.coordinates.sanitize_center(
            self.center, self.axis)
        xax = ds.coordinates.x_axis[self.axis]
        yax = ds.coordinates.y_axis[self.axis]
        bounds = (display_center[xax] - width[0] / 2,
                  display_center[xax] + width[0] / 2,
                  display_center[yax] - width[1] / 2,
                  display_center[yax] + width[1] / 2)

        return bounds

    def get_halo_map(self, ds, filename=None):
        """Get halo map and locate peaks

        Parameter
        ---------
        ds: yt.frontends.gadget.data_structures.GadgetHDF5Dataset

        Return
        ------
        peaklist: np.ndarray
            In which contains location of halo peaks.
        """
        # Get projected halo density
        halo_proj = ds.proj(self.halo_field, self.axis)
        # Get fiexed resolution buffer
        bounds = self.get_window_parameters(ds)
        halo_frb = FixedResolutionBuffer(halo_proj, bounds, self.imgsize)
        # YTarray to np.ndarray
        halo_map = np.array(halo_frb[self.halo_field], dtype=float)
        # Normalization
        halo_min = halo_map.min()
        halo_max = halo_map.max()
        halo_map_norm = (halo_map - halo_min) / (halo_max - halo_min)
        # Detect peaks
        # peaklist = self.locate_peak(halo_map_norm)

        if self.save == True:
            pz = yt.ProjectionPlot(ds, self.axis, self.halo_field,
                                   width=self.width)
            # pz.set_buff_size(self.buff_size)
            filepath = os.path.join(self.output_dir, filename)
            pz.save(filepath)

        return halo_map_norm

    def get_gas_map_draw(self, ds, filename=None):
        """Get 3D gas maps and mark peaks by self

        Parameter
        ---------
        ds = (ds_c1,ds_c2): yt.frontends.gadget.data_structures.
             GadgetHDF5Dataset

        peaklist: np.ndarray
            In which contains location of halo peaks

        """
        self.wlim = self.width[0]/2
        # Generate c1 gas map
        ds_c1 = ds[0]
        ad_c1 = ds_c1.all_data()
        # coordinates        
        x_c1 = np.array(ad_c1[('PartType0','particle_position_x')])
        y_c1 = np.array(ad_c1[('PartType0','particle_position_y')])
        z_c1 = np.array(ad_c1[('PartType0','particle_position_z')])
        den_c1 = np.array(ad_c1[('PartType0','Density')])
        # den_c1 = (den_c1 - den_c1.min())/(den_c1.max() - den_c1.min())
        # change units
        x_c1 = (x_c1 * self.unit_base['UnitLength_in_cm']) * 3.240779289469756e-25
        y_c1 = (y_c1 * self.unit_base['UnitLength_in_cm']) * 3.240779289469756e-25
        z_c1 = (z_c1 * self.unit_base['UnitLength_in_cm']) * 3.240779289469756e-25
        # Get indices
        idx_x = (x_c1 >= -self.wlim) * (x_c1 <= self.wlim) 
        idx_y = (y_c1 >= -self.wlim) * (y_c1 <= self.wlim)
        idx_z = (z_c1 >= -self.wlim) * (z_c1 <= self.wlim) 
        idx_c1 = idx_x * idx_y * idx_z
        # plot 
        # set axes and labels
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')
        ax.view_init(azim=-80,elev=20)
        ax.scatter3D(x_c1[idx_c1],y_c1[idx_c1],z_c1[idx_c1],s=0.01,
                    c= 'b', zdir='z',
                    marker='.',cmap=self.cmaps[0])
        # ax.axis([-box_lim,box_lim,-box_lim,box_lim,-box_lim,box_lim])
                
        # Generate c2 gas map
        ds_c2 = ds[1]
        ad_c2 = ds_c2.all_data()
        # coordinates        
        x_c2 = np.array(ad_c2[('PartType0','particle_position_x')])
        y_c2 = np.array(ad_c2[('PartType0','particle_position_y')])
        z_c2 = np.array(ad_c2[('PartType0','particle_position_z')])
        den_c2 = np.array(ad_c2[('PartType0','Density')])
        # den_c2 = (den_c2 - den_c2.min())/(den_c2.max() - den_c2.min())
        # change units
        x_c2 = (x_c2 * self.unit_base['UnitLength_in_cm']) * 3.240779289469756e-25
        y_c2 = (y_c2 * self.unit_base['UnitLength_in_cm']) * 3.240779289469756e-25
        z_c2 = (z_c2 * self.unit_base['UnitLength_in_cm']) * 3.240779289469756e-25
        # Get indices
        idx_x = (x_c2 >= -self.wlim) * (x_c2 <= self.wlim) 
        idx_y = (y_c2 >= -self.wlim) * (y_c2 <= self.wlim)
        idx_z = (z_c2 >= -self.wlim) * (z_c2 <= self.wlim) 
        idx_c2 = idx_x * idx_y * idx_z
        # plot 
        fig.hold()
        ax.scatter(x_c2[idx_c2],y_c2[idx_c2],z_c2[idx_c2],c= 'r',
                    zdir='z', marker='.',s=0.01)
                
        idx = re.findall(r'[0-9][0-9][0-9]', filename)
        i = int(idx[0])
        ax.text3D(-0.35, 0.35, 0.35, '%.2f Gyr' %
                         (i * 0.02))
        # mark
        mark_c1 = self.locate_peak(x_c1,y_c1,z_c1,den_c1)
        #print(mark_c1)
        ax.scatter3D(mark_c1[0],mark_c1[1],mark_c1[2], s=120,marker='x',c='r')
        mark_c2 = self.locate_peak(x_c2,y_c2,z_c2,den_c2)
        ax.scatter3D(mark_c2[0],mark_c2[1],mark_c2[2], s=80,marker='+',c='b')
        # set axes and labels
        box_lim = self.width[0]/2
        ax.set_xlabel('x (Mpc)')
        ax.set_ylabel('y (Mpc)')
        ax.set_zlabel('z (Mpc)')
        # set box limitation
        ax.set_xlim3d(-box_lim,box_lim)
        ax.set_ylim3d(-box_lim,box_lim)
        ax.set_zlim3d(-box_lim,box_lim)
        # plt.scatter(peaklist[1,1],peaklist[1,2],80,marker='x',c='r')
        # plt.scatter(peaklist[1,1],peaklist[1,2],80,marker='x',c='r')
        if self.save == True:
            filepath = os.path.join(self.output_dir, filename)
            fig.savefig(filepath)
            fig.clear()

    def locate_peak(self,x,y,z,density,width=2):
        """
        Locate peaks in the map
        """
        # Init
        wlim = width/2
        idx_x = (x >= -wlim) * (x <= wlim) 
        idx_y = (y >= -wlim) * (y <= wlim)
        idx_z = (z >= -wlim) * (z <= wlim) 
        indices = idx_x * idx_y * idx_z

        idx = np.where(density[indices]==density[indices].max())
        x = x[indices][idx]
        y = y[indices][idx]
        z = z[indices][idx]
        peaklist = [x,y,z]
        
        return peaklist

    def pix_to_unit(self, peaklist):
        """Locate the real pixels in the raw image."""
        rows, cols = self.imgsize
        pix_per_unit_col = self.width[0] / cols
        pix_per_unit_row = self.width[0] / rows
        peaklist[:, 1] = peaklist[:, 1] * pix_per_unit_col - self.width[0] / 2
        peaklist[:, 2] = peaklist[:, 2] * pix_per_unit_row - self.width[0] / 2

        return peaklist


def main(argv):
    """
    Prcessing on the whole simulation files.

    Parameter
    ---------
    argv: string
        argv = input_dir, output_dir
    """
    # Init
    # Input direction
    input_dir = argv[1]
    # Output direction
    output_dir = argv[2]
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    # Init GasTrack class
    gt = GasTrack(input_dir, output_dir)
    # Init sub particle hdft files
    ps_output = os.path.join(gt.input_dir, 'subs')
    if os.path.exists(ps_output):
        shutil.rmtree(ps_output)
        os.mkdir(ps_output)
    else:
        os.mkdir(ps_output)
    ps = partsplit.PartSplit(gt.input_dir, ps_output)
    # main circulation
    files = os.listdir(path=input_dir)
    files.sort()
    BeginIdx = 0
    EndinIdx = 500
    for i, f in enumerate(files):
        if os.path.splitext(f)[-1] == '.hdf5':
            snapidx = re.findall(r'[0-9][0-9][0-9]', f)
            snapidx = int(snapidx[0])
            if snapidx >= BeginIdx and snapidx <= EndinIdx:
                print('snap %03d' % snapidx)
                # Gen split particle hdf5 files
                h5name = ("snap_%03d.hdf5" % snapidx)
                ps.get_single_file(h5name)
                # Init save names
                haloname = ("halo_%03d_Projection_z_density_%f_%d_0.png" %
                                (snapidx, gt.width[0], gt.buff_size[0]))
                gasname = ("snap_%03d_Projection_z_density_%f_%d_0.png" %
                           (snapidx, gt.width[0], gt.buff_size[0]))
                # get peaks
                # gt.save = False
                width_old = gt.width
                gt.width = (2, 'Mpc')
                # cluster1
                f_c1 = ('subs/snap_%03d_c1.hdf5' % snapidx)
                ds1 = gt.load_hdf5(f_c1)
                # cluster2
                f_c2 = ('subs/snap_%03d_c2.hdf5' % snapidx)
                ds2 = gt.load_hdf5(f_c2)
                
                # Save results
                gt.width = width_old
                ds = gt.load_hdf5(f)

                gt.get_halo_map(ds, haloname)
                gt.get_gas_map_draw((ds1,ds2), gasname)
                os.remove(os.path.join(gt.input_dir,f_c1))
                os.remove(os.path.join(gt.input_dir,f_c2))
        else:
            pass

if __name__ == "__main__":
    # multi processing
    main(sys.argv)
