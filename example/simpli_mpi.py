import numpy as np
import h5py
import os

from mpi4py import MPI

import fastpli.simulation

FILE_NAME = os.path.abspath(__file__)
FILE_PATH = os.path.dirname(FILE_NAME)
FILE_BASE = os.path.basename(FILE_NAME)

with h5py.File('/tmp/fastpli.example.' + FILE_BASE + '.' +
               str(MPI.COMM_WORLD.Get_size()) + '.h5',
               'w',
               driver='mpio',
               comm=MPI.COMM_WORLD) as h5f:

    ### Setup Simpli for Tissue Generation
    simpli = fastpli.simulation.Simpli(MPI.COMM_WORLD)
    simpli.voxel_size = 1  # in mu meter
    simpli.set_voi([-60, 60, -60, 60, -30, 30])  # in mu meter
    simpli.ReadFiberFile(os.path.join(FILE_PATH, 'cube.h5'))
    simpli.fiber_bundles_properties = [[(0.333, -0.004, 10, 'p'),
                                        (0.666, 0, 5, 'b'),
                                        (1.0, 0.004, 1, 'r')]]

    print('VOI:', simpli.get_voi())
    print('Memory:', str(round(simpli.memory_usage('MB'), 2)) + ' MB')
    print(
        'Memory per thread:',
        str(round(simpli.memory_usage('MB') / MPI.COMM_WORLD.Get_size(), 2)) +
        ' MB')

    ### Generate Tissue
    print("Run Generation:")
    label_field, vec_field, tissue_properties = simpli.generate_tissue()

    simpli.save_mpi_array_as_h5(h5f, label_field.astype(np.uint16), 'tissue')
    simpli.save_mpi_array_as_h5(h5f, vec_field, 'vectorfield')

    ### Simulate PLI Measurement ###
    simpli.filter_rotations = np.deg2rad([0, 30, 60, 90, 120, 150])
    simpli.light_intensity = 26000  # a.u.
    simpli.untilt_sensor = True
    simpli.wavelength = 525  # in nm
    simpli.resolution = 20  # in mu meter
    TILTS = [(0, 0), (5.5, 0), (5.5, 90), (5.5, 180), (5.5, 270)]

    image_stack = [None] * len(TILTS)
    print("Run Simulation:")
    for t, (theta, phi) in enumerate(TILTS):
        image = simpli.run_simulation(label_field, vec_field, tissue_properties,
                                      np.deg2rad(theta), np.deg2rad(phi))

        simpli.save_mpi_array_as_h5(h5f, image, 'data/' + str(t), lock_dim=2)

# to apply optical resolution, the hole dataset has to be known
MPI.COMM_WORLD.barrier()

del label_field, vec_field  # free memory

if MPI.COMM_WORLD.Get_rank() == 0:
    with h5py.File(
            '/tmp/fastpli.example.' + FILE_BASE + '.' +
            str(MPI.COMM_WORLD.Get_size()) + '.h5', 'a') as h5f:

        with open(os.path.abspath(__file__), 'r') as f:
            h5f['script'] = f.read()

        print("Optic:")
        for t, (theta, phi) in enumerate(TILTS):

            # load data
            image = h5f['data/' + str(t)][:]

            # apply optic to simulation
            image = simpli.apply_optic(image, gain=3)
            h5f['optic/' + str(t)] = image

            # calculate modalities
            epa = simpli.apply_epa(image)
            h5f['epa/' + str(t) + '/transmittance'] = epa[0]
            h5f['epa/' + str(t) + '/direction'] = np.rad2deg(epa[1])
            h5f['epa/' + str(t) + '/retardation'] = epa[2]

            image_stack[t] = image

        print("Run ROFL analysis:")
        rofl_direction, rofl_incl, rofl_t_rel, _ = simpli.apply_rofl(
            image_stack, tilt_angle=np.deg2rad(5.5), gain=3)

        h5f['rofl/direction'] = np.rad2deg(rofl_direction)
        h5f['rofl/inclination'] = np.rad2deg(rofl_incl)
        h5f['rofl/trel'] = rofl_t_rel
