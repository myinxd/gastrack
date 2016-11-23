# !/usr/bin/python3
# Copyright (C) 2016 Zhixian MA <zxma_sjtu@qq.com>
"""
    A tool to count amont of particles in the provided regions.
"""

import os
import sys
import re
import h5py
import numpy as np

import utils


class Cluster:
    """
    The cluster class

    Parameters
    ----------
    filepath: string
        File path of the snap.
    clstype: string
        Major or minor cluster, can be 'maj' or 'min'.
    numhalo: int
        Number of halo particles in the major cluster.
    numgas: int
        Number of gas particles in the major cluster.

    Methods
    -------
    load_hdf5: load file
    get_idx: get indices of the cluster according to its type.
    get_particles: get particles in the cluster.
    """

    def __init__(self, filepath, clstype, numhalo, numgas):
        self.numhalo = numhalo
        self.numgas = numgas
        self.clstype = clstype
        self.filepath = filepath
        self.unit_len = 3.08568E21  # [cm]
        self.cm2Mpc = 3.240779289469756E-25
        # Open file
        try:
            self.load_hdf5()
        except IOError:
            print('File does not exist.')
            return

    def load_hdf5(self):
        """Open file"""
        self.filepart = h5py.File(self.filepath, 'r')

    def get_idx(self):
        """Get indices of the cluster"""
        # Init
        gas_group = self.filepart['PartType0']
        halo_group = self.filepart['PartType1']
        gas_id = np.array(gas_group['ParticleIDs'])
        halo_id = np.array(halo_group['ParticleIDs'])
        # judge cluster
        if self.clstype == 'maj':
            # major cluster
            self.gasidx = gas_id <= self.numgas
            self.haloidx = halo_id <= self.numhalo + len(gas_id)
            self.gas_id = gas_id[self.gasidx]
            self.halo_id = halo_id[self.haloidx]
        elif self.clstype == 'min':
            # minor cluster
            self.gasidx = gas_id > self.numgas
            self.haloidx = halo_id > self.numhalo + len(gas_id)
            self.gas_id = gas_id[self.gasidx]
            self.halo_id = halo_id[self.haloidx]

    def get_particles(self):
        """
        Get the cords and data

        Output
        ------
        cords: np.ndarray
            coordinates of the particles
        data: np.ndarray
            data with repect to the field
        """
        # Get indices
        self.get_idx()
        # Get cords and fields
        # gas
        gas_group = self.filepart['PartType0']
        # coordinates
        cords = np.array(gas_group['Coordinates'])
        self.gas_cords = cords[self.gasidx, :] * self.unit_len * self.cm2Mpc
        # field
        data = np.array(gas_group['Density'])
        self.gas_den = data[self.gasidx]
        # halo
        halo_group = self.filepart['PartType1']
        # coordinates
        cords = np.array(halo_group['Coordinates'])
        self.halo_cords = cords[self.haloidx, :] * self.unit_len * self.cm2Mpc
        # field
        data = np.array(halo_group['Masses'])
        self.halo_den = data[self.haloidx]


def get_cls(filepath, numhalo, numgas):
    """Get clusters

    Parameters
    ----------
    filepath: string
        File path of the snap.
    numhalo: int
        Number of halo particles in the major cluster.
    numgas: int
        Number of gas particles in the major cluster.

    Output
    ------
    maj_cls: Cluster object
        The major cluster.
    min_cls: Cluster object
        The minor cluster.
    """
    # Get clusters
    maj_cls = Cluster(filepath, 'maj', numhalo, numgas)
    min_cls = Cluster(filepath, 'min', numhalo, numgas)
    maj_cls.get_particles()
    min_cls.get_particles()

    return maj_cls, min_cls

def get_peak(maj_cls, step=0.005):
    """
    Generate maps, find peaks and save result.

    Parameters
    ----------
    maj_cls: Cluster object
        The major cluster
    step: double
        Step or width of the cubic for generating mosaic projection map.

    Output
    ------
    peak: list
        Coordinate of the peak in the projected map.
    """
    halo_z = utils.gen_mosaic(maj_cls.halo_cords, maj_cls.halo_den, step, 'z')
    peak_z = utils.get_peaks(halo_z, step)
    # y direction
    halo_y = utils.gen_mosaic(maj_cls.halo_cords, maj_cls.halo_den, step, 'y')
    peak_y = utils.get_peaks(halo_y, step)
    # x direction
    halo_x = utils.gen_mosaic(maj_cls.halo_cords, maj_cls.halo_den, step, 'x')
    peak_x = utils.get_peaks(halo_x, step)
    # Combine
    x = (peak_z[1] + peak_y[1]) / 2
    y = (peak_z[2] + peak_x[1]) / 2
    z = (peak_x[2] + peak_y[2]) / 2
    peak = [x, y, z]

    return peak

def calc_particles(maj_cls, min_cls, peak, reg_mode, reg_params):
    """
    Calculate amount of particles in the provided region

    Parameters
    ----------
    reg_mode: string
        Mode of the region, can be 'cir' or 'elp'
    reg_params: list
        Parameters of the region
    """
    # Init
    part_total = 0
    part_maj = 0
    part_min = 0

    # Calc particles
    if reg_mode == 'cir':
        # parmaters
        x_c = peak[0]
        y_c = peak[1]
        z_c = peak[2]
        radius = reg_params
        # maj
        maj_x = maj_cls.gas_cords[:, 0] - x_c
        maj_y = maj_cls.gas_cords[:, 1] - y_c
        maj_z = maj_cls.gas_cords[:, 2] - z_c
        # maj_z = maj_cls.gas_cords[:,2]
        maj_dist = np.sqrt(maj_x**2 + maj_y**2 + maj_z**2)
        maj_idx = maj_dist <= radius
        part_maj = maj_idx.sum()
        # min
        min_x = min_cls.gas_cords[:, 0] - x_c
        min_y = min_cls.gas_cords[:, 1] - y_c
        min_z = min_cls.gas_cords[:, 2] - z_c
        # min_z = min_cls.gas_cords[:,2]
        min_dist = np.sqrt(min_x**2 + min_y**2 + min_z**2)
        min_idx = min_dist <= radius
        part_min = min_idx.sum()
        # sum
        part_total = part_maj + part_min
    elif reg_mode == 'sec':
        # parmaters
        x_c = peak[0]
        y_c = peak[1]
        z_c = peak[2]
        radius_low = reg_params[0]
        radius_high = reg_params[1]
        angle_low = reg_params[2]
        angle_high = reg_params[3]

        if angle_low >= 2 * np.pi:
            angle_low -= 2 * np.pi
            angle_high -= 2 * np.pi

        if angle_high > 2 * np.pi:
            angle_low = [angle_low, 0]
            angle_high = [2 * np.pi, angle_high - 2 * np.pi]
        else:
            angle_low = [angle_low]
            angle_high = [angle_high]

        maj_idx = np.zeros(maj_cls.gas_id.shape)
        maj_idx = maj_idx.astype(bool)
        min_idx = np.zeros(min_cls.gas_id.shape)
        min_idx = min_idx.astype(bool)

        for i in range(len(angle_low)):
            # maj
            maj_x = maj_cls.gas_cords[:, 0] - x_c
            maj_y = maj_cls.gas_cords[:, 1] - y_c
            maj_z = maj_cls.gas_cords[:, 2] - z_c
            # maj_z = maj_cls.gas_cords[:,2]
            maj_dist = np.sqrt(maj_x**2 + maj_y**2 + maj_z**2)
            maj_ang = np.arcsin(np.abs(maj_y) / maj_dist)
            maj_idx_dist = (maj_dist <= radius_high) * (maj_dist >= radius_low)
            # Quarant1
            idx_q1 = (maj_x >= 0) * (maj_y >= 0)
            idx_ang1 = (maj_ang <= angle_high[i]) * (maj_ang >= angle_low[i])
            # Quarant2
            idx_q2 = (maj_x < 0) * (maj_y >= 0)
            idx_ang2 = ((np.pi - maj_ang) <=
                        angle_high[i]) * ((np.pi - maj_ang) >= angle_low[i])
            # Quarant3
            idx_q3 = (maj_x < 0) * (maj_y < 0)
            idx_ang3 = ((np.pi + maj_ang) <=
                        angle_high[i]) * ((np.pi + maj_ang) >= angle_low[i])
            # Quarant4
            idx_q4 = (maj_x >= 0) * (maj_y <= 0)
            idx_ang4 = ((2 * np.pi - maj_ang) <=
                        angle_high[i]) * ((2 * np.pi - maj_ang) >= angle_low[i])
            # Combine idx
            maj_idx_t = (idx_ang1 * idx_q1) + (idx_ang2 * idx_q2) + \
                (idx_ang3 * idx_q3) + (idx_ang4 * idx_q4)
            maj_idx_t = (maj_idx_dist) * (maj_idx_t)
            maj_idx = maj_idx + maj_idx_t
            part_maj += maj_idx_t.sum()
            # min
            min_x = min_cls.gas_cords[:, 0] - x_c
            min_y = min_cls.gas_cords[:, 1] - y_c
            min_z = min_cls.gas_cords[:, 2] - z_c
            min_dist = np.sqrt(min_x**2 + min_y**2 + min_z**2)
            min_ang = np.arcsin(np.abs(min_y) / min_dist)
            min_idx_dist = (min_dist <= radius_high) * (min_dist >= radius_low)
            # Quarant1
            idx_q1 = (min_x >= 0) * (min_y >= 0)
            idx_ang1 = (min_ang <= angle_high[i]) * (min_ang >= angle_low[i])
            # Quarant2
            idx_q2 = (min_x < 0) * (min_y >= 0)
            idx_ang2 = ((np.pi - min_ang) <=
                        angle_high[i]) * ((np.pi - min_ang) >= angle_low[i])
            # Quarant3
            idx_q3 = (min_x < 0) * (min_y < 0)
            idx_ang3 = ((np.pi + min_ang) <=
                        angle_high[i]) * ((np.pi + min_ang) >= angle_low[i])
            # Quarant4
            idx_q4 = (min_x >= 0) * (min_y <= 0)
            idx_ang4 = ((2 * np.pi - min_ang) <=
                        angle_high[i]) * ((2 * np.pi - min_ang) >= angle_low[i])
            # Combine idx
            min_idx_t = (idx_ang1 * idx_q1) + (idx_ang2 * idx_q2) + \
                (idx_ang3 * idx_q3) + (idx_ang4 * idx_q4)
            min_idx_t = (min_idx_dist) * (min_idx_t)
            min_idx = min_idx + min_idx_t
            part_min += min_idx_t.sum()
        # sum
        part_total = part_maj + part_min

    else:
        print("Mode %s is not supported at present" % reg_mode)

    partlist = [part_total, part_maj, part_min]

    return partlist, maj_idx, min_idx

def main(argv):
    """The main method"""
    # Init
    file1 = argv[1]
    file2 = argv[2]
    # get id
    snapid1 = re.findall(r'[0-9][0-9][0-9]', file1)
    snapid1 = int(snapid1[0])
    snapid2 = re.findall(r'[0-9][0-9][0-9]', file2)
    snapid2 = int(snapid2[0])

    outpath = argv[3]
    # Init of parameters
    numhalo = int(argv[4])  # 734866
    numgas = int(argv[5])  # 704860
    step = 0.01

    # get cls of file1
    maj_cls_f1, min_cls_f1 = get_cls(file1, numhalo, numgas)
    # get cls of file2
    maj_cls_f2, min_cls_f2 = get_cls(file2, numhalo, numgas)
    # get peak
    peak1 = get_peak(maj_cls_f1, step)
    peak2 = get_peak(maj_cls_f2, step)
    # Calc particles of file2
    print('Searching for particles at %.2f Gyr ...' % (snapid2 * 0.02))
    part_ori, maj_idx_f2, min_idx_f2 = calc_particles(maj_cls_f2,
                                                      min_cls_f2,
                                                      peak2, 'cir', 50 / 1000)
    print('Total particles at %.2f Gyr: %d' % (snapid2 * 0.02, part_ori[0]))
    print('Major particles at %.2f Gyr: %d' % (snapid2 * 0.02, part_ori[1]))
    print('Minor particles at %.2f Gyr: %d' % (snapid2 * 0.02, part_ori[2]))

    # Calc particles of file1
    # calc circles
    print('Searching for particles in the circular region at %.2f Gyr...'
          % (snapid1 * 0.02))
    reg_mode = 'sec'
    reg_params = [float(argv[6]), float(argv[7]),
                  float(argv[8]) / 180 * np.pi,
                  float(argv[9]) / 180 * np.pi]
    part_sec, maj_idx_c, min_idx_c = calc_particles(
        maj_cls_f1, min_cls_f1, peak1, reg_mode, reg_params)
    print('Total particles at %.2f Gyr: %d' % (snapid1 * 0.02, part_sec[0]))
    print('Major particles at %.2f Gyr: %d' % (snapid1 * 0.02, part_sec[1]))
    print('Minor particles at %.2f Gyr: %d' % (snapid1 * 0.02, part_sec[2]))
    # diff
    part_maj_c = utils.cmp_id(maj_cls_f1, maj_cls_f2, maj_idx_c, maj_idx_f2)
    print('Major particles from %.2f Gyr: %d' % (snapid2 * 0.02, part_maj_c))
    part_min_c = utils.cmp_id(min_cls_f1, min_cls_f2, min_idx_c, min_idx_f2)
    print('Minor particles from %.2f Gyr: %d' % (snapid2 * 0.02, part_min_c))

    # save
    # filename = os.path.join(outdir,'particles.txt')
    filename = outpath
    if os.path.exists(filename):
        os.remove(filename)
        f = open(filename, 'a')
    else:
        f = open(filename, 'a')

    f.write("Particles at %0.2f Gyr.\n" % (snapid2 * 0.02))
    f.write('Total particles at %.2f Gyr: %d\n' %
            (snapid2 * 0.02, part_ori[0]))
    f.write('Major particles at %.2f Gyr: %d\n' %
            (snapid2 * 0.02, part_ori[1]))
    f.write('Minor particles at %.2f Gyr: %d\n' %
            (snapid2 * 0.02, part_ori[2]))
    f.write('\n')

    f.write("Particles at %0.2f Gyr in the section region.\n" %
            (snapid1 * 0.02))
    f.write('Total particles at %.2f Gyr: %d\n' %
            (snapid1 * 0.02, part_sec[0]))
    f.write('Major particles at %.2f Gyr: %d\n' %
            (snapid1 * 0.02, part_sec[1]))
    f.write('Minor particles at %.2f Gyr: %d\n' %
            (snapid1 * 0.02, part_sec[2]))
    f.write('Major particles from %.2f Gyr: %d\n' %
            (snapid2 * 0.02, part_maj_c))
    f.write('Minor particles from %.2f Gyr: %d\n' %
            (snapid2 * 0.02, part_min_c))

    f.close()

if __name__ == "__main__":
    main(sys.argv)
