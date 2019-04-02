import fastpli.io
import fastpli.model
import fastpli.simulation

import h5py
import numpy as np
import os

INPUT_FILE = '/localdata/data/xy-crossing_90deg.dat'
OUTPUT_FILE_FIBER = '/tmp/test.dat'
OUTPUT_FILE_TISSUE = '/tmp/test.h5'
# FILE_PATH = os.path.dirname(os.path.abspath(__file__))

# read fiber_bundle file
print('loading files')
fiber_bundles = fastpli.io.fiber.read(INPUT_FILE)

# rnd movement
for fiber_bundle in fiber_bundles:
    for fiber in fiber_bundle:
        fiber += np.random.uniform(0, 0.1, fiber.shape)

# setup solver
print('setting up solver')
solver = fastpli.model.Solver()
solver.fiber_bundles = fiber_bundles
solver.drag = 0
solver.obj_min_radius = 5
solver.obj_mean_length = 1
solver.omp_num_threads = 8

# run solver
print('running solver')
i = 0
MAX_SOVLER_STEPS = 200

for i in range(MAX_SOVLER_STEPS):
    solved = solver.step()

    if solved:
        break

    if i % 10 == 0:
        print("step:", i, solver.num_obj, solver.num_col_obj)
        solver.draw_scene()

# save data
if OUTPUT_FILE_FIBER:
    fastpli.io.fiber.save(OUTPUT_FILE_FIBER, solver.fiber_bundles)

# tissue generation
simpli = fastpli.simulation.Simpli()
simpli.debug = False

simpli.pixel_size = 0.125
simpli.dim = [400, 400, 240]
simpli.dim_origin = -simpli.dim * simpli.pixel_size / 2.0  # position of coordinate (0,0,0) in mu-meter
print("simpli.voi:", simpli.voi)
# simpli.voi in mu-meter, can also be setted, dim and origin will be changed accordingly
# e.g. simpli.voi = [25, 75, 25, 75, 15, 45]

simpli.fiber_bundles = solver.fiber_bundles
simpli.fiber_bundles_properties = [[(1, 0, 0, 'p')]]

label_field, _, _ = simpli.GenerateTissue(only_label=True)

if OUTPUT_FILE_TISSUE:
    with h5py.File(OUTPUT_FILE_TISSUE, 'w') as h5f:
        dset = h5f.create_dataset(
            'tissue', simpli.dim, dtype=np.uint16, compression="gzip")
        dset[:] = label_field