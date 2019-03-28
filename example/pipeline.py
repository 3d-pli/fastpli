import fastpli.io
import fastpli.model

import numpy as np
import time
import os

FILE_PATH = os.path.dirname(os.path.abspath(__file__))

# read fiber_bundle file
print('loading files')
fiber_bundles = fastpli.io.fiber.read('/localdata/data/xy-crossing_90deg.dat')

# setup solver
print('setting up solver')
solver = fastpli.model.Solver()
solver.fiber_bundles = fiber_bundles
solver.drag = 0
solver.obj_min_radius = 0
solver.obj_mean_length = 15
solver.omp_num_threads = 8

# run solver
print('running solver')
i = 0
MAX_SOVLER_STEPS = 100
while not solver.step() and i < MAX_SOVLER_STEPS:

    if i % 10 == 0:
        print("step:", i, solver.num_obj, solver.num_col_obj)
        solver.draw_scene()

    i += 1
