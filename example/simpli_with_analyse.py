import h5py
import numpy as np
import sys
import os

import scipy.misc
import nibabel as nib

from fastpli.simulation import Simpli
from fastpli.analysis import images

FILE_PATH = os.path.dirname(os.path.abspath(__file__))
FILE_NAME = os.path.splitext(os.path.basename(__file__))[0]
FILE_OUTPUT = '/tmp/'

# FILE_PATH = os.path.dirname(os.path.abspath(__file__))

# Read fiber data and prepair for PliGenerator
# TODO: write json -> parameter function

simpli = Simpli()
# PliGeneration ###
simpli.pixel_size = 4
simpli.dim = [160, 160, 15]
simpli.fiber_bundles = [[np.array([[0, 0, 30, 120], [1000, 1000, 30, 120]])]]
simpli.fiber_bundles_properties = [[(1.0, 0.004, 1, 'p')]]

print("MemoryUseage:", '~' + str(int(simpli.MemoryUseage())) + ' MB')

with h5py.File(os.path.join(FILE_OUTPUT, 'example.' + FILE_NAME + '.h5'),
               'w') as h5f:

    with open(os.path.abspath(__file__), 'r') as f:
        h5f['script'] = f.read()

    label_field, vec_field, tissue_properties = simpli.GenerateTissue()

    dset = h5f.create_dataset('tissue', label_field.shape, np.uint16)
    dset[:] = label_field
    dset.attrs['element_size_um'] = simpli.pixel_size

    h5f['vectorfield'] = vec_field
    h5f['vectorfield'].attrs['element_size_um'] = simpli.pixel_size

    # PliSimulation ###
    simpli.filter_rotations = np.deg2rad([0, 30, 60, 90, 120, 150])
    simpli.light_intensity = 26000
    simpli.pixel_size = 1
    simpli.untilt_sensor = True
    simpli.wavelength = 525

    tilts = [(0, 0), (5.5, 0), (5.5, 90), (5.5, 180), (5.5, 270)]
    image_stack = []

    print("Run Simulations:")
    for t, (theta, phi) in enumerate(tilts):
        print('Step:', t, theta, phi)
        image = simpli.RunSimulation(label_field, vec_field, tissue_properties,
                                     np.deg2rad(theta), np.deg2rad(phi))
        h5f['data/' + str(t)] = image
        h5f['data/' + str(t)].attrs['element_size_um'] = simpli.pixel_size
        h5f['data/' + str(t)].attrs['theta'] = theta
        h5f['data/' + str(t)].attrs['phi'] = phi

        opt_pixel_size = 64.0
        image = simpli.apply_optic(image, opt_pixel_size)
        h5f['optic/' + str(t)] = image
        h5f['optic/' + str(t)].attrs['element_size_um'] = opt_pixel_size
        h5f['optic/' + str(t)].attrs['theta'] = theta
        h5f['optic/' + str(t)].attrs['phi'] = phi

        image_stack.append(image)

    print("Run ROFL:")
    rofl_direction, rofl_incl, rofl_t_rel = simpli.apply_rofl(image_stack)
    h5f['rofl/direction'] = np.rad2deg(rofl_direction)
    h5f['rofl/inclination'] = np.rad2deg(rofl_incl)
    h5f['rofl/t_rel'] = rofl_t_rel

    h5f['rofl/direction'].attrs['element_size_um'] = opt_pixel_size
    h5f['rofl/inclination'].attrs['element_size_um'] = opt_pixel_size
    h5f['rofl/t_rel'].attrs['element_size_um'] = opt_pixel_size

    print("Unit Vectors")
    unit_x, unit_y, unit_z = images.unit_vectors(rofl_direction, rofl_incl)
    img = nib.Nifti1Image(unit_x, np.eye(4))
    nib.save(img,
             os.path.join(FILE_OUTPUT, 'example.' + FILE_NAME + '.UnitX.nii'))
    img = nib.Nifti1Image(unit_y, np.eye(4))
    nib.save(img,
             os.path.join(FILE_OUTPUT, 'example.' + FILE_NAME + '.UnitY.nii'))
    img = nib.Nifti1Image(unit_z, np.eye(4))
    nib.save(img,
             os.path.join(FILE_OUTPUT, 'example.' + FILE_NAME + '.UnitZ.nii'))

    print("FOMs")
    img = images.fom_hsv_black(rofl_direction, rofl_incl)
    scipy.misc.imsave(
        os.path.join(FILE_OUTPUT, 'example.' + FILE_NAME + '.hsv_black.png'),
        np.swapaxes(img, 0, 1))

    img = images.fom_rgb(rofl_direction, rofl_incl)
    scipy.misc.imsave(
        os.path.join(FILE_OUTPUT, 'example.' + FILE_NAME + '.rgb.png'),
        np.swapaxes(img, 0, 1))

    img = images.hsvblack_sphere()
    print(img.shape)
    print(type(img))
    print(type(img.flatten()[0]))
    scipy.misc.imsave(
        os.path.join(FILE_OUTPUT,
                     'example.' + FILE_NAME + '.hsv_black_sphere.png'),
        np.swapaxes(img, 0, 1))

    img = images.rgb_sphere()
    scipy.misc.imsave(
        os.path.join(FILE_OUTPUT, 'example.' + FILE_NAME + '.rgb_sphere.png'),
        np.swapaxes(img, 0, 1))
