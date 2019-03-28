import fastpli.io
import fastpli.model

import numpy as np
import os

FILE_PATH = os.path.dirname(os.path.abspath(__file__))

# read fiber_bundle file
print('loading files')
fiber_bundles = fastpli.io.fiber.read('/localdata/data/xy-crossing_90deg.dat')

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
solver.obj_mean_length = 2
solver.omp_num_threads = 8

# run solver
print('running solver')
i = 0
MAX_SOVLER_STEPS = 2

for i in range(MAX_SOVLER_STEPS):
    solved = solver.step()

    if solved:
        break

    if i % 10 == 0:
        print("step:", i, solver.num_obj, solver.num_col_obj)
        solver.draw_scene()

# save data
fastpli.io.fiber.save('/tmp/test.dat', solver.fiber_bundles)
