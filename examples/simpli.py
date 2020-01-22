import fastpli.simulation
import fastpli.analysis
import fastpli.tools
import fastpli.io

import numpy as np
import h5py
import os

import imageio

np.random.seed(42)

FILE_NAME = os.path.abspath(__file__)
FILE_PATH = os.path.dirname(FILE_NAME)
FILE_BASE = os.path.basename(FILE_NAME)

with h5py.File('/tmp/fastpli.example.' + FILE_BASE + '.h5', 'w') as h5f:
    # save script
    h5f['version'] = fastpli.__version__
    with open(os.path.abspath(__file__), 'r') as f:
        h5f['parameter/script'] = f.read()
        h5f['parameter/pip_freeze'] = fastpli.tools.helper.pip_freeze()

    # Setup Simpli for Tissue Generation
    simpli = fastpli.simulation.Simpli()
    simpli.omp_num_threads = 2
    simpli.voxel_size = 0.5  # in mu meter
    simpli.set_voi([-60] * 3, [60] * 3)  # in mu meter
    simpli.fiber_bundles = fastpli.io.fiber.load(
        os.path.join(FILE_PATH, 'cube.dat'))

    # (layer_scale, dn, mu, optical-axis)
    # b: background, p: parallel, r: radial
    simpli.fiber_bundles_properties = [[(0.333, -0.004, 10, 'p'),
                                        (0.666, 0, 5, 'b'),
                                        (1.0, 0.004, 1, 'r')]]

    print('VOI:', simpli.get_voi())
    print('Memory:', str(round(simpli.memory_usage('MB'), 2)) + ' MB')

    # Generate Tissue
    print("Run Generation:")
    label_field, vec_field, tissue_properties = simpli.generate_tissue()

    h5f['tissue/label_field'] = label_field.astype(np.uint16)
    h5f['tissue/vectorfield'] = vec_field
    h5f['tissue/tissue_properties'] = tissue_properties

    # Simulate PLI Measurement
    simpli.filter_rotations = np.deg2rad([0, 30, 60, 90, 120, 150])
    simpli.light_intensity = 26000  # a.u.
    simpli.interpolate = True
    simpli.wavelength = 525  # in nm
    simpli.resolution = 20  # in mu meter
    simpli.sensor_gain = 3
    simpli.optical_sigma = 0.71  # in voxel size
    simpli.tilts = np.deg2rad([(0, 0), (5.5, 0), (5.5, 90), (5.5, 180),
                               (5.5, 270)])

    tilting_stack = [None] * 5
    print("Run Simulation:")
    for t, (theta, phi) in enumerate(simpli.tilts):
        print(round(np.rad2deg(theta), 1), round(np.rad2deg(phi), 1))
        images = simpli.run_simulation(label_field, vec_field,
                                       tissue_properties, theta, phi)

        h5f['simulation/data/' + str(t)] = images

        # apply optic to simulation
        print(images.shape)
        images = simpli.apply_optic(images)
        h5f['simulation/optic/' + str(t)] = images

        # calculate modalities
        epa = simpli.apply_epa(images)
        h5f['analysis/epa/' + str(t) + '/transmittance'] = epa[0]
        h5f['analysis/epa/' + str(t) + '/direction'] = np.rad2deg(epa[1])
        h5f['analysis/epa/' + str(t) + '/retardation'] = epa[2]

        tilting_stack[t] = images

    # save mask for analysis
    mask = np.sum(label_field, 2) > 0
    mask = simpli.apply_optic_reshape(1.0 * mask) > 0.1
    h5f['simulation/optic/mask'] = np.uint8(mask)
    mask = None  # keep analysing all pixels

    print("Run ROFL analysis:")
    rofl_direction, rofl_incl, rofl_t_rel, _ = simpli.apply_rofl(tilting_stack,
                                                                 mask=mask)

    h5f['analysis/rofl/direction'] = np.rad2deg(rofl_direction)
    h5f['analysis/rofl/inclination'] = np.rad2deg(rofl_incl)
    h5f['analysis/rofl/trel'] = rofl_t_rel

    # def data2image(data):
    #     return np.swapaxes(np.flip(data, 1), 0, 1)

    # imageio.imwrite(
    #     os.path.join(FILE_PATH, 'simpli.png'),
    #     data2image(
    #         fastpli.analysis.images.fom_hsv_black(rofl_direction, rofl_incl)))
