# mpirun -n $i python3 -m mpi4py examples/simpli_mpi.py

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

with h5py.File(f'{FILE_OUT}_{MPI.COMM_WORLD.Get_size()}.h5',
               'w',
               driver='mpio',
               comm=MPI.COMM_WORLD) as h5f:

    # saveing string not supported in mpi yet
    # h5f['version'] = fastpli.__version__
    # h5f['parameter/pip_freeze'] = fastpli.tools.helper.pip_freeze()
    # with open(os.path.abspath(__file__), 'r') as f:
    #     h5f['parameter/script'] = f.read()

    # Setup Simpli for Tissue Generation
    simpli = fastpli.simulation.Simpli(MPI.COMM_WORLD)
    simpli.omp_num_threads = 1
    simpli.voxel_size = 0.5  # in micro meter
    # simpli.set_voi([-20, -5, -10], [20, 5, 10])  # in micro meter
    simpli.set_voi([-60] * 3, [60] * 3)  # in micro meter
    simpli.fiber_bundles = fastpli.io.fiber_bundles.load(
        os.path.join(FILE_PATH, 'cube.dat'))

    # define layers (e.g. axon, myelin) inside fibers of each fiber_bundle fiber_bundle
    simpli.fiber_bundles_properties = [[(0.333, -0.004, 10, 'p'),
                                        (0.666, 0, 5, 'b'),
                                        (1.0, 0.004, 1, 'r')]]
    # (_0, _1, _2, _3)
    # _0: layer_scale times radius
    # _1: strength of birefringence
    # _2: absorption coefficient mu: I = I*exp(-mu*x)
    # _3: model: 'p'-parallel, 'r'-radial or 'b'-background

    print('VOI:', simpli.get_voi())
    print('Memory:', str(round(simpli.memory_usage('MB'), 2)) + ' MB')

    # Generate Tissue
    print('Run Generation:')
    tissue, optical_axis, tissue_properties = simpli.generate_tissue()

    simpli.mpi.save_h5(h5f, 'tissue/tissue', tissue.astype(np.uint16))
    simpli.mpi.save_h5(h5f, 'tissue/optical_axis', optical_axis)
    h5f['tissue/tissue_properties'] = tissue_properties

    # Simulate PLI Measurement
    simpli.filter_rotations = np.deg2rad([0, 30, 60, 90, 120, 150])
    simpli.light_intensity = 26000  # a.u.
    simpli.interpolate = "Slerp"
    # simpli.untilt_sensor_view = True
    simpli.wavelength = 525  # in nm
    simpli.pixel_size = 20  # in micro meter
    simpli.optical_sigma = 0.71  # in pixel size
    simpli.noise_model = lambda x: np.random.negative_binomial(
        x / (3 - 1), 1 / 3)
    simpli.tilts = np.deg2rad([(0, 0), (5.5, 0), (5.5, 90), (5.5, 180),
                               (5.5, 270)])

    tilting_stack = [None] * 5
    print('Run Simulation:')
    for t, (theta, phi) in enumerate(simpli.tilts):
        print(f'theta: {np.rad2deg(theta):.1f}, phi: {np.rad2deg(phi):.1f}')
        images = simpli.run_simulation(tissue, optical_axis, tissue_properties,
                                       theta, phi)
        simpli.mpi.save_image_h5(h5f, f'simulation/data/{t}', images)

        # apply optic to simulation
        if MPI.COMM_WORLD.Get_rank() == 0:
            images = simpli.apply_optic(images)  #use mp_pool to speed up
            dim = np.array(images.shape)
        else:
            dim = np.empty(3, dtype=np.int64)
        MPI.COMM_WORLD.Bcast(dim, root=0)
        if MPI.COMM_WORLD.Get_rank() > 0:
            images = np.empty(dim, np.float64)
        MPI.COMM_WORLD.Bcast(images, root=0)

        images = simpli.mpi.split_optic(images)  # splitting tasks
        simpli.mpi.save_split_h5(h5f, f'simulation/optic/{t}', images)

        # calculate modalities
        epa = simpli.apply_epa(images)
        simpli.mpi.save_split_h5(h5f, f'analysis/epa/{t}/transmittance', epa[0])
        simpli.mpi.save_split_h5(h5f, f'analysis/epa/{t}/direction',
                                 np.rad2deg(epa[1]))
        simpli.mpi.save_split_h5(h5f, f'analysis/epa/{t}/retardation', epa[2])

        tilting_stack[t] = images.copy()

    print('Run ROFL analysis:')
    rofl_direction, rofl_incl, rofl_t_rel, _ = simpli.apply_rofl(tilting_stack)

    simpli.mpi.save_split_h5(h5f, f'analysis/rofl/direction',
                             np.rad2deg(rofl_direction))
    simpli.mpi.save_split_h5(h5f, f'analysis/rofl/inclination',
                             np.rad2deg(rofl_incl))
    simpli.mpi.save_split_h5(h5f, f'analysis/rofl/trel', rofl_t_rel)

    print('Done')
    print('You can look at the data e.g with Fiji and the hdf5 plugin')
