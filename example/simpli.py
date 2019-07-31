import numpy as np
import h5py
import os

import fastpli.simulation

FILE_NAME = os.path.abspath(__file__)
FILE_PATH = os.path.dirname(FILE_NAME)
FILE_BASE = os.path.basename(FILE_NAME)

with h5py.File('/tmp/fastpli.example.' + FILE_BASE + '.h5', 'w') as h5f:
    ### save script ###
    with open(os.path.abspath(__file__), 'r') as f:
        h5f['script'] = f.read()

    ### Setup Simpli for Tissue Generation
    simpli = fastpli.simulation.Simpli()
    simpli.pixel_size = 1
    simpli.dim = [100, 100, 100]
    simpli.ReadFiberFile(os.path.join(FILE_PATH, 'cube.h5'))
    simpli.fiber_bundles_properties = [[(0.333, -0.004, 10, 'p'),
                                        (0.666, 0, 5, 'b'),
                                        (1.0, 0.004, 1, 'r')]]

    print('VOI:', simpli.get_voi())
    print('Memory:', str(round(simpli.MemoryUseage('MB'), 2)) + ' MB')

    ### Generate Tissue
    label_field, vec_field, tissue_properties = simpli.GenerateTissue()

    h5f['tissue'] = label_field.astype(np.uint16)
    h5f['vectorfield'] = vec_field

    ### Simulate PLI Measurement ###
    simpli.filter_rotations = np.deg2rad([0, 30, 60, 90, 120, 150])
    simpli.light_intensity = 26000  # a.u.
    simpli.pixel_size = 1  # mu-m
    simpli.untilt_sensor = True
    simpli.wavelength = 525  # nm
    TILTS = [(0, 0), (5.5, 0), (5.5, 90), (5.5, 180), (5.5, 270)]

    for t, (theta, phi) in enumerate(TILTS):
        print(theta, phi)
        image = simpli.RunSimulation(label_field, vec_field, tissue_properties,
                                     np.deg2rad(theta), np.deg2rad(phi))
        h5f['data/' + str(t)] = image
