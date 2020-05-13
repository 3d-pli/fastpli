import fastpli.simulation
import fastpli.analysis
import fastpli.tools
import fastpli.io

from mpi4py import MPI
import numpy as np
import h5py
import os

import imageio

np.random.seed(42)

FILE_NAME = os.path.abspath(__file__)
FILE_PATH = os.path.dirname(FILE_NAME)
FILE_BASE = os.path.basename(FILE_NAME)
FILE_OUT = os.path.join(FILE_PATH, f'fastpli.example.{FILE_BASE}')

print(f'creating file: {FILE_OUT}_{MPI.COMM_WORLD.Get_size()}.h5')

# Setup Simpli for Tissue Generation

# with h5py.File(f'{FILE_OUT}_{MPI.COMM_WORLD.Get_size()}.h5',
#                'w',
#                driver='mpio',
#                comm=MPI.COMM_WORLD) as h5f:
with h5py.File(
        f'{FILE_OUT}_{MPI.COMM_WORLD.Get_size()}_{MPI.COMM_WORLD.Get_rank()}.h5',
        'w') as h5f:
    simpli = fastpli.simulation.Simpli(MPI.COMM_WORLD)
    simpli.omp_num_threads = 1
    simpli.voxel_size = 0.5  # in µm meter
    # simpli.set_voi([-20, -5, -10], [20, 5, 10])  # in µm meter
    simpli.set_voi([-60] * 3, [60] * 3)  # in µm meter
    simpli.fiber_bundles = fastpli.io.fiber_bundles.load(
        os.path.join(FILE_PATH, 'cube.dat'))

    # define layers (e.g. axon, myelin) inside fibers of each fiber_bundle fiber_bundle
    simpli.fiber_bundles_properties = [[(0.333, -0.004, 10, 'p'),
                                        (0.666, 0, 5, 'b'),
                                        (1.0, 0.004, 1, 'r')]]
    # (_0, _1, _2, _3)
    # _0: layer_scale times radius
    # _1: strength of birefringence
    # _2: absorption coefficient µ: I = I*exp(-µ*x)
    # _3: model: 'p'-parallel, 'r'-radial or 'b'-background

    print('VOI:', simpli.get_voi())
    print('Memory:', str(round(simpli.memory_usage('MB'), 2)) + ' MB')

    # Generate Tissue
    print('Run Generation:')
    tissue, optical_axis, tissue_properties = simpli.generate_tissue()

    h5f['tissue/tissue'] = tissue.astype(np.uint16)
    h5f['tissue/optical_axis'] = optical_axis

    # Simulate PLI Measurement
    simpli.filter_rotations = np.deg2rad([0, 30, 60, 90, 120, 150])
    simpli.light_intensity = 26000  # a.u.
    simpli.interpolate = True
    # simpli.untilt_sensor_view = True
    simpli.wavelength = 525  # in nm
    simpli.pixel_size = 20  # in µm meter
    simpli.sensor_gain = 3
    simpli.optical_sigma = 0.71  # in voxel size
    simpli.tilts = np.deg2rad([(0, 0), (5.5, 0), (5.5, 90), (5.5, 180),
                               (5.5, 270)])

    tilting_stack = [None] * 5
    print('Run Simulation:')
    for t, (theta, phi) in enumerate(simpli.tilts):
        print(f'theta: {np.rad2deg(theta):.1f}, phi: {np.rad2deg(phi):.1f}')
        images = simpli.run_simulation(tissue, optical_axis, tissue_properties,
                                       theta, phi)

        h5f['simulation/data/' + str(t)] = images

# to apply optical resolution, the hole dataset has to be known
MPI.COMM_WORLD.barrier()
del optical_axis  # free memory

# TRANSMITE IMAGES to every rank?

# if MPI.COMM_WORLD.Get_rank() == 0:
#     with h5py.File(f'{FILE_OUT}_{MPI.COMM_WORLD.Get_size()}.h5', 'a') as h5f:

#         # save script
#         h5f['version'] = fastpli.__version__
#         h5f['parameter/pip_freeze'] = fastpli.tools.helper.pip_freeze()
#         h5f['tissue/tissue_properties'] = tissue_properties
#         with open(os.path.abspath(__file__), 'r') as f:
#             h5f['parameter/script'] = f.read()

#         print("Run Optic:")
#         for t, (theta, phi) in enumerate(simpli.tilts):

#             # load data
#             images = h5f['simulation/data/' + str(t)][:]

#             # apply optic to simulation
#             images = simpli.apply_optic(images)
#             h5f['simulation/optic/' + str(t)] = images

#             # calculate modalities
#             epa = simpli.apply_epa(images)
#             h5f['analysis/epa/' + str(t) + '/transmittance'] = epa[0]
#             h5f['analysis/epa/' + str(t) + '/direction'] = np.rad2deg(epa[1])
#             h5f['analysis/epa/' + str(t) + '/retardation'] = epa[2]

#             tilting_stack[t] = images

#         # save mask for analysis
#         mask = np.sum(tissue, 2) > 0
#         mask = simpli.apply_optic_resample(1.0 * mask) > 0.1
#         h5f['simulation/optic/mask'] = np.uint8(mask)
#         mask = None  # keep analysing all pixels

#         # TODO: This can be parallel again
#         print('Run ROFL analysis:')
#         rofl_direction, rofl_incl, rofl_t_rel, _ = simpli.apply_rofl(
#             tilting_stack, mask=mask)

#         h5f['analysis/rofl/direction'] = np.rad2deg(rofl_direction)
#         h5f['analysis/rofl/inclination'] = np.rad2deg(rofl_incl)
#         h5f['analysis/rofl/trel'] = rofl_t_rel

#         def data2image(data):
#             return np.swapaxes(np.flip(data, 1), 0, 1)

#         print(f'creating Fiber Orientation Map: {FILE_OUT}.png')

#         imageio.imwrite(
#             f'{FILE_OUT}.png',
#             data2image(
#                 fastpli.analysis.images.fom_hsv_black(rofl_direction,
#                                                       rofl_incl)))

#     print('Done')
#     print('You can look at the data e.g with Fiji and the hdf5 plugin')
