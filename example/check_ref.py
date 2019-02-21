import os
import h5py
import numpy as np
from fastpli.simulation import Simpli

FILE_PATH = os.path.dirname(os.path.abspath(__file__))

# MPI.Init()
simpli = Simpli()
# PliGeneration ###
simpli.pixel_size = 1
simpli.dim = [100, 100, 100]
simpli.dim_origin = [0, 0, 0]

simpli.ReadFiberFile('example/cube.h5')
simpli.SetFiberProperties([[(0.333, 0.004, 10, 'p'), (
    0.666, -0.004, 5, 'b'), (1.0, 0.004, 1, 'r')]])

# manipulation of fibers
simpli.RotateVolumeAroundPoint(np.deg2rad(
    20), np.deg2rad(-10), np.deg2rad(5), [10, -5, 7.5])
simpli.TranslateVolume([25, -15, 50])

with h5py.File(os.path.join(FILE_PATH, 'ref_new.h5'), 'w') as h5f:

    label_field, _, tissue_properties = simpli.GenerateTissue(only_label=True)
    label_field, vec_field, tissue_properties = simpli.GenerateTissue()

    h5f['tissue'] = label_field
    h5f['vectorfield'] = vec_field

    # PliSimulation ###
    simpli.setup.filter_rotations = np.deg2rad([0, 30, 60, 90, 120, 150])
    simpli.setup.light_intensity = 26000
    simpli.setup.resolution = 1
    simpli.setup.untilt_sensor = True
    simpli.setup.wavelength = 525

    # simpli.InitSimulation(label_field, vec_field, tissue_properties)

    print("run_simulation: 0")
    h5f['data/0'] = simpli.run_simulation(
        label_field,
        vec_field,
     tissue_properties,
     0,
     0)

    print("run_simulation: 1")
    h5f['data/1'] = simpli.run_simulation(
        label_field,
        vec_field,
        tissue_properties,
        np.deg2rad(5.5),
        np.deg2rad(0))

    print("run_simulation: 2")
    h5f['data/2'] = simpli.run_simulation(
        label_field,
        vec_field,
        tissue_properties,
        np.deg2rad(5.5),
        np.deg2rad(90))

    print("run_simulation: 3")
    h5f['data/3'] = simpli.run_simulation(
        label_field,
        vec_field,
        tissue_properties,
        np.deg2rad(5.5),
        np.deg2rad(180))

    print("run_simulation: 4")
    h5f['data/4'] = simpli.run_simulation(
        label_field,
        vec_field,
        tissue_properties,
        np.deg2rad(5.5),
        np.deg2rad(270))


os.system(
    'h5diff ' + os.path.join(
        FILE_PATH,
         'ref_new.h5') + ' ' + os.path.join(
             FILE_PATH,
             'ref_68ae085.h5'))
