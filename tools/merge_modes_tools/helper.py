#!/bin/python3
import glob
import os
import h5py
import numpy as np
import sys
import getopt


def get_max_mode(path = './'):

    real_dirs = glob.glob(path + 'Re*')
    imag_dirs = glob.glob(path + 'Im*')
    num_real_dirs = 0
    num_imag_dirs = 0

    for obj in real_dirs:
        if os.path.isdir(obj):
            num_real_dirs += 1

    for obj in imag_dirs:
        if os.path.isdir(obj):
            num_imag_dirs += 1

    assert num_real_dirs - 1 == num_imag_dirs, 'Wrong number of folders!'

    return num_real_dirs - 1


def get_file_list(path = './'):

    file_list = glob.glob(path + 'Re0/*_*.h5')
    file_list = [os.path.basename(file) for file in file_list]
    return file_list
    # e.g.: file_list = ['charge_00000000.h5', 'charge_00000001.h5']

def merge(path = './',phi_degree=0):
    phi = float(phi_degree) / 180 * np.pi
    max_mode = get_max_mode(path)
    file_list = get_file_list(path)
    num_files = len(file_list)
    if not os.path.exists(path+'Merged_angle_'+str(phi_degree)):
        os.mkdir(path+'Merged_angle_'+str(phi_degree))

    filename = path + 'Re0/' + file_list[0]
    with h5py.File(filename, 'r') as h5file:

        # read dataset name and attributes
        dset_name = list(h5file.keys())[1]
        h5dset = h5file[dset_name]
        dset_attrs = {name: value for name, value in h5dset.attrs.items()}
        dset_size = h5dset.shape

        # read axes and their attributes
        h5axis1 = h5file['AXIS']['AXIS1']
        h5axis2 = h5file['AXIS']['AXIS2']
        axis1 = np.array(h5axis1)
        axis2 = np.array(h5axis2)
        axis1[0] = -axis1[1]
        axis1_attrs = {name: value for name, value in h5axis1.attrs.items()}
        axis2_attrs = {name: value for name, value in h5axis2.attrs.items()}

    for filename in file_list:

        print('Merging file ' + path + filename + '...')

        # initialize merged dataset
        merged_dset = np.zeros((dset_size[0], dset_size[1] * 2 - 1))

        # add m = 0 mode to the merged dataset
        with h5py.File(path + 'Re0/' + filename, 'r') as h5file:

            # read file attributes
            file_attrs = {name: value for name, value in h5file.attrs.items()}

            dset = np.array(h5file[dset_name])
            merged_dset[:, 0:dset_size[1] - 1] += np.fliplr(dset[:, 1:])
            merged_dset[:, dset_size[1] - 1:] += dset

        # add m > 0 modes to the merged dataset
        for mode in np.arange(1, max_mode + 1):

            with h5py.File(path + 'Re' + str(mode) + '/' + filename, 'r') as h5file:
                dset = np.array(h5file[dset_name])
                merged_dset[:, 0:dset_size[1] - 1] += np.fliplr(dset[:, 1:]) * 2 * np.cos(mode * (phi + np.pi))
                merged_dset[:, dset_size[1] - 1:] += dset * 2 * np.cos(mode * phi)

            with h5py.File(path + 'Im' + str(mode) + '/' + filename, 'r') as h5file:
                dset = np.array(h5file[dset_name])
                merged_dset[:, 0:dset_size[1] - 1] -= np.fliplr(dset[:, 1:]) * np.sin(mode * (phi + np.pi))
                merged_dset[:, dset_size[1] - 1:] -= dset * np.sin(mode * phi)

        # write to the h5 files
        with h5py.File(path + 'Merged_angle_'+str(phi_degree)+'/' + filename, 'w') as h5file:

            # write file attributes
            for key, val in file_attrs.items():
                h5file.attrs[key] = val

            # create new dataset and write the attributes
            h5dset = h5file.create_dataset(dset_name, data=merged_dset, dtype='f8')
            for key, val in dset_attrs.items():
                h5dset.attrs[key] = val

            # create group for axis and write the informations
            h5axis = h5file.create_group('AXIS')
            h5axis1 = h5axis.create_dataset('AXIS1', data=axis1)
            h5axis2 = h5axis.create_dataset('AXIS2', data=axis2)
            for key, val in axis1_attrs.items():
                h5axis1.attrs[key] = val
            for key, val in axis2_attrs.items():
                h5axis2.attrs[key] = val
