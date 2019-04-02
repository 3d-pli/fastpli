import fastpli.io
import fastpli.model
import fastpli.simulation

import h5py
import numpy as np
import os

np.random.seed(42)

INPUT_FILE = '/localdata/data/xy-crossing_90deg.dat'
OUTPUT_FILE_FIBER = '/tmp/test_chunked.dat'
OUTPUT_FILE_TISSUE = '/tmp/test_chunked.h5'
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
MAX_SOVLER_STEPS = 10

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
# simpli.dim = [400, 400, 240]
# simpli.dim_origin = -simpli.dim * simpli.pixel_size / 2.0  # position of coordinate (0,0,0) in mu-meter
print("simpli.voi:", simpli.voi)
# simpli.voi in mu-meter, can also be setted, dim and origin will be changed accordingly
# e.g. simpli.voi = [25, 75, 25, 75, 15, 45]

simpli.voi = [25, 75, 25, 75, 15, 45]

simpli.fiber_bundles = solver.fiber_bundles
simpli.fiber_bundles_properties = [[(1, 0, 0, 'p')]]

xvoi = (simpli.voi[0], simpli.voi[1])
X_STEP = 5

if xvoi[0] % X_STEP != 0 or xvoi[1] % X_STEP != 0:
    raise ValueError("not implemented yet")

print(simpli.dim)
print(simpli.dim_origin)
print(simpli.voi)

with h5py.File(OUTPUT_FILE_TISSUE, 'w') as h5f:
    dset = h5f.create_dataset(
        'tissue', simpli.dim, dtype=np.uint16, compression="gzip")

    x_list = [(x, x + X_STEP) for x in range(xvoi[0], xvoi[1], X_STEP)]

    # print(x_list)

    for x0, x1 in x_list:
        simpli.voi = [x0, x1, 25, 75, 15, 45]

        xdim = (int((x0 - xvoi[0]) // simpli.pixel_size),
                int((x1 - xvoi[0]) // simpli.pixel_size))

        label_field, _, _ = simpli.GenerateTissue(only_label=True)

        print(simpli.dim)
        print(simpli.dim_origin)
        print(simpli.voi)

        # print(simpli.dim)
        # print(simpli.voi)
        # print(x0, x1)
        # print(xdim[0], xdim[1])
        # print(label_field.shape)

        dset[xdim[0]:xdim[1], :, :] = label_field
