#!/usr/bin/python3
# Author: Zhixian MA <zxma_sjtu@qq.com>
#
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
"""

import os
import re
import sys
import shutil
import numpy as np
import yt
from yt.visualization.fixed_resolution import FixedResolutionBuffer

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
                 buffsize=(3200, 3200)):
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.save = True
        # Projection axis
        self.axis = 'z'
        # center
        self.center = 'c'
        # width
        self.width = width
        # buff_size
        self.buff_size = buffsize
        # image size
        self.imgsize = (800, 800)
        # unit_base
        self.unit_base = {'UnitLength_in_cm': 3.08568E21,
                          'UnitMass_in_g': 1.989E43,
                          'UnitVelocity_in_cm_per_s': 100000}
        # Bound box
        bbox_lim = 1E5  # [kpc]
        self.bbox = [[-bbox_lim, bbox_lim],
                     [-bbox_lim, bbox_lim],
                     [-bbox_lim, bbox_lim]]
        # fields
        self.halo_field = ('deposit', 'PartType1_density')
        self.gas_field = ('gas', 'density')

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
            pz.set_buff_size(self.buff_size)
            filepath = os.path.join(self.output_dir, filename)
            pz.save(filepath)

        return halo_map_norm

    def get_gas_map(self, ds, peaklist, filename=None):
        """Get gas map and mark peaks.

        Parameter
        ---------
        ds: yt.frontends.gadget.data_structures.GadgetHDF5Dataset

        peaklist: np.ndarray
            In which contains location of halo peaks

        """
        # Generate gas map
        pz = yt.ProjectionPlot(ds, self.axis, self.gas_field,
                               width=self.width)
        # Set buff_size
        pz.set_buff_size(self.buff_size)
        # Markers
        if peaklist.shape[0] == 1:
            # one peak
            pz.annotate_marker((peaklist[0, 1], peaklist[0, 2]),
                               coord_system='plot',
                               plot_args={'color': 'blue', 's': 500})
        else:
            pz.annotate_marker((peaklist[0, 1], peaklist[0, 2]),
                               coord_system='plot',
                               plot_args={'color': 'blue', 's': 500})
            pz.annotate_marker((peaklist[1, 1], peaklist[1, 2]),
                               marker="+",
                               coord_system='plot',
                               plot_args={'color': 'red', 's': 300})
        idx = re.findall(r'[0-9][0-9][0-9]', filename)
        i = int(idx[0])
        pz.annotate_text((-0.3, 0.3), '%.2f Gyr' %
                         (i * 0.02), coord_system='plot')
        if self.save == True:
            filepath = os.path.join(self.output_dir, filename)
            pz.save(filepath)

    def locate_peak(self, img_mat):
        """
        Locate peaks in the map

        References
        ----------
        [1] http://stackoverflow.com/questions/3684484/peak-detection-
        in-a-2d-array
        [2] http://stackoverflow.com/questions/9111711/get-coordinates-
        of-local-maxima-in-2d-array-above-certain-value
        """
        # Init
        peaks = []
        cord_x = []
        cord_y = []
        rows, cols = img_mat.shape
        # Find peaks
        peak_max = img_mat.max()
        peak_y, peak_x = np.where(img_mat == peak_max)
        peak_y = int(round(np.mean(peak_y)))
        peak_x = int(round(np.mean(peak_x)))
        # Judege and fill
        # append
        peaks.append(peak_max)
        cord_x.append(peak_x)
        cord_y.append(peak_y)

        peaklist = np.array([peaks, cord_x, cord_y]).transpose()
        print(peaklist)
        # pix to unit
        peaklist = self.pix_to_unit(peaklist)

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
    BeginIdx = 3
    for i, f in enumerate(files):
        if os.path.splitext(f)[-1] == '.hdf5':
            snapidx = re.findall(r'[0-9][0-9][0-9]', f)
            snapidx = int(snapidx[0])
            if snapidx == BeginIdx:
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
                gt.save = True
                width_old = gt.width
                gt.width = (2, 'Mpc')
                # cluster1
                f_c1 = ('subs/snap_%03d_c1.hdf5' % snapidx)
                ds1 = gt.load_hdf5(f_c1)
                halo_c1 = ("halo_%03d_Projection_z_density_%f_%d_c1.png" %
                               (snapidx, gt.width[0], gt.buff_size[0]))
                halo_map_c1 = gt.get_halo_map(ds1, halo_c1)
                peak1 = gt.locate_peak(halo_map_c1)
                os.remove(os.path.join(gt.input_dir, f_c1))
                # cluster2
                f_c2 = ('subs/snap_%03d_c2.hdf5' % snapidx)
                ds2 = gt.load_hdf5(f_c2)
                halo_c2 = ("halo_%03d_Projection_z_density_%f_%d_c2.png" %
                           (snapidx, gt.width[0], gt.buff_size[0]))
                halo_map_c2 = gt.get_halo_map(ds2, halo_c2)
                peak2 = gt.locate_peak(halo_map_c2)
                os.remove(os.path.join(gt.input_dir, f_c2))
                # Combine
                peaklist = np.row_stack((peak1, peak2))
                # get all images
                gt.save = True
                gt.width = width_old
                ds = gt.load_hdf5(f)
                halo_map = gt.get_halo_map(ds, haloname)
                gt.get_gas_map(ds, peaklist, gasname)
        else:
            pass

if __name__ == "__main__":
    main(sys.argv)
