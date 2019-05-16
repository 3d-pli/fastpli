import fastpli.model
import fastpli.simulation
import fastpli.io
import fastpli.analysis

import scipy.misc
import nibabel as nib
import h5py

import numpy as np
import time
import os

import apsg

FILE_PATH = os.path.dirname(os.path.abspath(__file__))
FILE_NAME = os.path.splitext(os.path.basename(__file__))[0]
FILE_OUTPUT = '/tmp/'

np.random.seed(42)
''' create fiber_bundle(s)
'''
# create fiber_bundle trajectory
fiber_bundles = [[]]
# fiber_bundles[-1].append([[0, 0, 0, 1], [10, 10, 10, 1]])
# fiber_bundles[-1].append([[0, 0, 0, 1], [0, 10, 10, 1]])

kent = apsg.helpers.KentDistribution((1, 0, 0), (0, 1, 0), (0, 0, 1), 16, 5)

# v = np.random.uniform(-100, 100, (100, 3))
v = kent.rvs(1000) * 300
for i in range(v.shape[0]):
    offset = np.random.uniform(-100, 100, (3,))
    fiber_bundles[-1].append(
        [[-v[i, 0] + offset[0], -v[i, 1] + offset[1], -v[i, 2] + offset[2], 5],
         [v[i, 0] + offset[0], v[i, 1] + offset[1], v[i, 2] + offset[2], 5]])

# print(fiber_bundles)

fastpli.io.fiber.save('/tmp/fastpli.example.model_solver.dat', fiber_bundles)
''' setup solver
'''
solver = fastpli.model.Solver()
solver.fiber_bundles = fiber_bundles
solver.drag = 0
solver.obj_min_radius = 0
solver.obj_mean_length = 15
solver.omp_num_threads = 8

solved_fbs = solver.fiber_bundles
fastpli.io.fiber.save('/tmp/fastpli.example.model_solver_.dat', solved_fbs)
''' run solver
'''
for i in range(250):

    solved = solver.step()

    if i % 25 == 0:
        print("step:", i, solver.num_obj, solver.num_col_obj)
        solver.draw_scene()

    if solved:
        print("step:", i, solver.num_obj, solver.num_col_obj)
        solver.draw_scene()
        break

# PliGeneration ###
simpli = fastpli.simulation.Simpli()
simpli.omp_threads(2)
simpli.pixel_size = 1
simpli.resolution = 20
simpli.set_voi([-250, 250, -250, 250, -30, 30])

simpli.fiber_bundles = fiber_bundles
simpli.fiber_bundles_properties = [[(1.0, -0.001, 1, 'p')]
                                  ]  # scale, dn, mu, orientation

print("MemoryUseage:", '~' + str(int(simpli.MemoryUseage())) + ' MB')

with h5py.File(os.path.join(FILE_OUTPUT, 'example.' + FILE_NAME + '.new.h5'),
               'w') as h5f:

    with open(os.path.abspath(__file__), 'r') as f:
        h5f['script'] = f.read()

    label_field, vec_field, tissue_properties = simpli.GenerateTissue()

    dset = h5f.create_dataset('tissue', label_field.shape, np.uint16)
    dset[:] = label_field

    h5f['vectorfield'] = vec_field

    # PliSimulation ###
    simpli.filter_rotations = np.deg2rad([0, 30, 60, 90, 120, 150])
    simpli.light_intensity = 26000
    simpli.untilt_sensor = True
    simpli.wavelength = 525

    tilts = [(0, 0), (5.5, 0), (5.5, 90), (5.5, 180), (5.5, 270)]

    image_stack = []

    mask = None
    # mask = np.sum(label_field, 2) > 0
    # mask = simpli.apply_resize_mask(mask) > 0.1

    print("Run Simulations:")
    for t, (theta, phi) in enumerate(tilts):
        print('Step:', t, theta, phi)
        image = simpli.RunSimulation(label_field, vec_field, tissue_properties,
                                     np.deg2rad(theta), np.deg2rad(phi))
        h5f['data/' + str(t)] = image

        image = simpli.apply_optic(image, gain=3)
        h5f['optic/' + str(t)] = image

        epa = simpli.apply_epa(image)
        h5f['epa/' + str(t) + '/transmittance'] = epa[0]
        h5f['epa/' + str(t) + '/direction'] = np.rad2deg(epa[1])
        h5f['epa/' + str(t) + '/retardation'] = epa[2]

        image_stack.append(image)

    # save mask for analysis, but restore None value so that everything gets analysed
    mask = np.sum(label_field, 2) > 0
    mask = simpli.apply_resize_mask(mask) > 0.1
    h5f['optic/mask'] = np.uint8(mask)
    mask = None  # analyse all pixels

    print("Run ROFL:")
    rofl_direction, rofl_incl, rofl_t_rel, _ = simpli.apply_rofl(
        image_stack, tilt_angle=np.deg2rad(5.5), gain=3, mask=mask)

    h5f['rofl/direction'] = np.rad2deg(rofl_direction)
    h5f['rofl/inclination'] = np.rad2deg(rofl_incl)
    h5f['rofl/trel'] = rofl_t_rel

    print("Unit Vectors")
    unit_x, unit_y, unit_z = fastpli.analysis.images.unit_vectors(
        rofl_direction, rofl_incl, mask)
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
    img = fastpli.analysis.images.fom_hsv_black(rofl_direction, rofl_incl, mask)
    scipy.misc.imsave(
        os.path.join(FILE_OUTPUT, 'example.' + FILE_NAME + '.hsv_black.png'),
        np.swapaxes(np.flip(img, 1), 0, 1))

    img = fastpli.analysis.images.fom_rgb(rofl_direction, rofl_incl, mask)
    scipy.misc.imsave(
        os.path.join(FILE_OUTPUT, 'example.' + FILE_NAME + '.rgb.png'),
        np.swapaxes(np.flip(img, 1), 0, 1))

    img = fastpli.analysis.images.hsvblack_sphere()
    scipy.misc.imsave(
        os.path.join(FILE_OUTPUT,
                     'example.' + FILE_NAME + '.hsv_black_sphere.png'),
        np.swapaxes(np.flip(img, 1), 0, 1))

    img = fastpli.analysis.images.rgb_sphere()
    scipy.misc.imsave(
        os.path.join(FILE_OUTPUT, 'example.' + FILE_NAME + '.rgb_sphere.png'),
        np.swapaxes(np.flip(img, 1), 0, 1))

    test = np.zeros((256, 256, 3), np.uint8)
    test[0, 0, :] = [0, 0, 255]
    test[-1, 0, :] = [255, 0, 0]
    test[0, -1, :] = [0, 255, 0]
    test[-1, -1, :] = [255, 255, 0]
    scipy.misc.imsave(
        os.path.join(FILE_OUTPUT, 'example.' + FILE_NAME + '.test.png'),
        np.swapaxes(np.flip(test, 1), 0, 1))
