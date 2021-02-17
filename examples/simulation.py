import fastpli.simulation
import fastpli.io

import multiprocessing as mp
import matplotlib.pyplot as plt
import numpy as np
import h5py
import os

pool = mp.Pool(2)
np.random.seed(42)

FILE_NAME = os.path.abspath(__file__)
FILE_PATH = os.path.dirname(FILE_NAME)

# Setup Simpli for Tissue Generation
simpli = fastpli.simulation.Simpli()
simpli.omp_num_threads = 2

# define volume
simpli.voxel_size = 5  # in micro meter
simpli.set_voi([-2000, -2000, -50], [2000, 2000, 50])  # in micro meter
print('Memory:', str(round(simpli.memory_usage('MB'), 2)) + ' MB')

# define model
simpli.fiber_bundles = fastpli.io.fiber_bundles.load(
    os.path.join(FILE_PATH, 'optic_chiasm.dat'))
simpli.fiber_bundles.layers = [[(1.0, -0.0001, 5, 'p')]] * len(
    simpli.fiber_bundles)
# (_0, _1, _2, _3)
# _0: layer_scale times radius
# _1: strength of birefringence
# _2: absorption coefficient I = I*exp(-mu*x)
# _3: model: 'p'-parallel, 'r'-radial or 'b'-background

# define pli setup
simpli.filter_rotations = np.deg2rad([0, 30, 60, 90, 120, 150])
simpli.light_intensity = 26000  # a.u.
simpli.wavelength = 525  # in nm
simpli.pixel_size = 20  # in micro meter
simpli.optical_sigma = 0.71  # in pixel size
simpli.noise_model = lambda x: np.random.negative_binomial(x / (3 - 1), 1 / 3)
simpli.tilts = np.deg2rad([(0, 0), (5.5, 0), (5.5, 90), (5.5, 180), (5.5, 270)])

print(f'creating file: simpli.h5')
with h5py.File(f'simpli.h5', 'w') as h5f:
    with open(os.path.abspath(__file__), 'r') as script:
        results = simpli.run_pipeline(h5f=h5f,
                                      script=script.read(),
                                      save=['all'],
                                      crop_tilt=True,
                                      mp_pool=pool)

        tissue, optical_axis, tissue_properties = results[0]
        tilting_stack = results[1]
        rofl_direction, rofl_incl, rofl_t_rel = results[2]
        fom = results[3]

    fig, axs = plt.subplots(1, 4, figsize=(15, 5))
    axs[0].imshow(np.swapaxes(rofl_direction, 0, 1))
    axs[1].imshow(np.swapaxes(rofl_incl, 0, 1))
    axs[2].imshow(np.swapaxes(rofl_t_rel, 0, 1))
    axs[3].imshow(np.swapaxes(fom, 0, 1))
    plt.show()

pool.close()
