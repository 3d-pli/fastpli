import time
import h5py
import numpy as np
import sys
import os

from fastpli.simulation import Simpli

FILE_PATH = os.path.dirname(os.path.abspath(__file__))


# Read fiber data and prepair for PliGenerator
# TODO: write json -> parameter function

global_start = time.time()

simpli = Simpli()
# PliGeneration ###
simpli.pixel_size = 1
simpli.dim = [100, 100, 100]
simpli.ReadFiberFile('example/cube.h5')
simpli.SetFiberProperties([[(0.333, 0.004, 10, 'p'), (
    0.666, -0.004, 5, 'b'), (1.0, 0.004, 1, 'r')]])

# manipulation of fibers
simpli.RotateVolumeAroundPoint(np.deg2rad(
    20), np.deg2rad(-10), np.deg2rad(5), [10, -5, 7.5])
simpli.TranslateVolume([25, -15, 50])

with h5py.File(os.path.join(FILE_PATH, 'output.h5'), 'w') as h5f:

    start = time.time()
    label_field, _, tissue_properties = simpli.GenerateTissue(only_label=True)
    label_field, vec_field, tissue_properties = simpli.GenerateTissue()

    h5f['tissue'] = label_field
    h5f['vectorfield'] = vec_field
    end = time.time()

    print("TissueGeneration:", end - start)

    # PliSimulation ###
    simpli.setup.filter_rotations = np.deg2rad([0, 30, 60, 90, 120, 150])
    simpli.setup.light_intensity = 26000
    simpli.setup.resolution = 1
    simpli.setup.untilt_sensor = True
    simpli.setup.wavelength = 525

    start = time.time()
    print("run_simulation: 0")
    image = simpli.run_simulation(
        label_field, vec_field, tissue_properties, 0, 0)
    h5f['data/0'] = image

    print("run_simulation: 1")
    image = simpli.run_simulation(
        label_field,
        vec_field,
        tissue_properties,
        np.deg2rad(5.5),
        np.deg2rad(0))
    h5f['data/1'] = image

    print("run_simulation: 2")
    image = simpli.run_simulation(
        label_field,
        vec_field,
        tissue_properties,
        np.deg2rad(5.5),
        np.deg2rad(90))
    h5f['data/2'] = image

    print("run_simulation: 3")
    image = simpli.run_simulation(
        label_field,
        vec_field,
        tissue_properties,
        np.deg2rad(5.5),
        np.deg2rad(180))
    h5f['data/3'] = image

    print("run_simulation: 4")
    image = simpli.run_simulation(
        label_field,
        vec_field,
        tissue_properties,
        np.deg2rad(5.5),
        np.deg2rad(270))
    h5f['data/4'] = image

    end = time.time()
    print("run_simulation:", end - start)

global_end = time.time()
print("GlobalRuntime:", global_end - global_start)
