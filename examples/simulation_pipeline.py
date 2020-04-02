import numpy as np
import copy
import h5py
import os
import sys
import glob

import multiprocessing as mp
pool = mp.Pool(2)

np.random.seed(42)

import fastpli.simulation
import fastpli.io

# reproducability

FILE_NAME = os.path.abspath(__file__)
FILE_PATH = os.path.dirname(FILE_NAME)
FILE_BASE = os.path.basename(FILE_NAME)

# Setup Simpli for Tissue Generation
simpli = fastpli.simulation.Simpli()
simpli.omp_num_threads = 2

# define model
simpli.voxel_size = 0.5  # in mu meter
simpli.set_voi([-60] * 3, [60] * 3)  # in mu meter
simpli.fiber_bundles = fastpli.io.fiber_bundles.load(
    os.path.join(FILE_PATH, 'cube.dat'))

# define layers (e.g. axon, myelin) inside fibers of each fiber_bundle fiber_bundle
simpli.fiber_bundles_properties = [[(0.333, -0.004, 10, 'p'),
                                    (0.666, 0, 5, 'b'), (1.0, 0.004, 1, 'r')]]
# (_0, _1, _2, _3)
# _0: layer_scale times radius
# _1: strength of birefringence
# _2: absorption coefficient I = I*exp(-mu*x)
# _3: model: 'p'-parallel, 'r'-radial or 'b'-background

# define pli setup
simpli.filter_rotations = np.deg2rad([0, 30, 60, 90, 120, 150])  # in deg
simpli.light_intensity = 26000  # a.u.
simpli.interpolate = True
simpli.wavelength = 525  # in nm
simpli.resolution = 5  # in mu meter
simpli.tilts = np.deg2rad(
    np.array([(0, 0), (5.5, 0), (5.5, 90), (5.5, 180), (5.5, 270)]))  # in deg
simpli.sensor_gain = 3
simpli.optical_sigma = 0.71  # in voxel size
simpli.verbose = True

file_name = 'fastpli.example.' + FILE_BASE + '.h5'
print(f"creating file: {file_name}")

with h5py.File(file_name, 'w') as h5f:
    with open(os.path.abspath(__file__), 'r') as script:
        simpli.run_pipeline(h5f=h5f,
                            script=script.read(),
                            save=["all"],
                            crop_tilt=True,
                            mp_pool=pool)

pool.close()

print("Done - You can look at the data e.g with Fiji and the hdf5 plugin")
