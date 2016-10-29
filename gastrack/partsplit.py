# Author Zhixian MA <zxma_sjtu@qq.com>

"""
This module split simulated galaxy particles hdf5 file into two sub
hdf5 files.
"""

import os
import re
import h5py
import numpy as np


class PartSplit:
    """
    The partticle split class

    Parameters
    ----------
    input_dir: string
        The direction holding those father hdf5 files.
    output_dir: string
        The direction to save output files.
    numpart_main: int
        Number of particles of the main cluster.

    Methods
    -------
    load_hdf5
        load the father hdf5 file.
    get_sub_part:
        Get subset of particles according to the masses ratio.

    References
    ----------
    [1] Collette, A.,
        "Python and HDF5",
        O'reilly, 2013.
    [2] h5py doc
        https://docs.h5py.org
    """

    def __init__(self, input_dir, output_dir, numpart_main=66985):
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.numpart_main = numpart_main

    def load_hdf5(self, fname):
        """
        Open the hdf5 file with h5py and return the object.

        Paramater
        ---------
        fname: string
            File name

        Return
        ------
        filepart: hdf5 object
            The object holding particle datasets.
        """
        # open
        filepath = os.path.join(self.input_dir, fname)
        filepart = h5py.File(filepath, 'r')

        return filepart

    def get_sub_part(self, filepart, partidx, fname):
        """
        Split the whole particles into two subsets according to
        the number of particles for each.

        Parameters
        ----------
        filepart: hdf5 object
            The object holding particle detasets.
        partidx: np.ndarray list
            Indices of the halo and gas particles to be extracted.
        fname: string
            File names of the sub set.
        """
        # Init
        filepath = os.path.join(self.output_dir, fname)
        if os.path.exists(filepath):
            os.remove(filepath)
        subpart = h5py.File(filepath, 'a')
        # Header
        header = subpart.create_group('Header')
        for key, value in filepart['Header'].attrs.items():
            header.attrs[key] = value
        # Change NumPart_Thisfile
        numpart = partidx[0].shape[0]
        # header.attrs['NumPaat_Total'] = [numpart,numpart,0,0,0,0]
        header.attrs['NumPart_ThisFile'] = [numpart, numpart, 0, 0, 0, 0]

        # PartType0
        parttype0 = subpart.create_group('PartType0')
        group_type0 = filepart['PartType0']
        for key in group_type0.keys():
            dataset = group_type0[key]
            # Extract
            partset = dataset[:]
            try:
                partset = partset[partidx[0], :]
            except IndexError:
                partset = partset[partidx[0]]
            parttype0.create_dataset(key, data=partset)

        # PartType1
        parttype1 = subpart.create_group('PartType1')
        group_type1 = filepart['PartType1']
        for key in group_type1.keys():
            dataset = group_type1[key]
            # Extract
            partset = dataset[:]
            try:
                partset = partset[partidx[1], :]
            except IndexError:
                partset = partset[partidx[1]]
            parttype1.create_dataset(key, data=partset)

        # close
        subpart.close()

    def get_single_file(self, fname):
        """
        Processing on single file.
        """
        # Init
        # output
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
        # Circle
        snapidx = re.findall(r'[0-9][0-9][0-9]', fname)
        snapidx = int(snapidx[0])
        print('Snap %03d' % snapidx)
        # load file
        filepart = self.load_hdf5(fname)
        # cluster1
        # Indices
        # gas
        particle_gas = filepart['/PartType0/ParticleIDs']
        gas_idx = np.where(particle_gas[:] < self.numpart_main)[0]
        # halo
        particle_halo = filepart['/PartType1/ParticleIDs']
        numgas = particle_gas[:].shape[0]
        halo_idx = np.where(particle_halo[:] <
                            self.numpart_main + numgas)[0]
        # idx
        partidx_c1 = (gas_idx, halo_idx)
        # split
        subname_c1 = ("snap_%03d_c1.hdf5" % snapidx)
        self.get_sub_part(filepart, partidx_c1, subname_c1)

        # cluster2
        # Indices
        # gas
        gas_idx = np.where(particle_gas[:] > self.numpart_main)[0]
        # halo
        halo_idx = np.where(particle_halo[:] >
                            self.numpart_main + numgas)[0]
        # idx
        partidx_c2 = (gas_idx, halo_idx)
        # split
        subname_c2 = ("snap_%03d_c2.hdf5" % snapidx)
        self.get_sub_part(filepart, partidx_c2, subname_c2)

        filepart.close()

        def get_multi_files(self):
            """
            Processing on multiple files.
            """
            # Init
            files = os.listdir(self.input_dir)
            files.sort()
            # output
            if not os.path.exists(self.output_dir):
                os.mkdir(self.output_dir)
            # Circle
            for f in files:
                if os.path.splitext(f)[-1] == '.hdf5':
                    self.get_single_file(f)
