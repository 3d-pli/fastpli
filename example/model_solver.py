import fastpli
from fastpli.model import bundle_generator
import fastpli.io

import numpy as np
import time
import os

FILE_PATH = os.path.dirname(os.path.abspath(__file__))
np.random.seed(42)
''' create fiber_bundle(s)
'''
# create fiber_bundle trajectory
fiber_bundle_trj = []
for phi in range(0, 91, 5):
    fiber_bundle_trj.append(
        [200 * np.sin(np.deg2rad(phi)), 200 * np.cos(np.deg2rad(phi)), 0, 1])
fiber_bundles = [[]]
fiber_bundles[-1].append(np.array(fiber_bundle_trj))

# fill fiber_bundle with fibers from a 2d population
filled_fiber_bundles = [[]]
population = bundle_generator.populate_circle(64, 1.5 * 5)
for fb in fiber_bundles:
    for fiber_obj in fb:
        fibers = bundle_generator.populate_object(fiber_obj, population, 5)
        for f in fibers:
            filled_fiber_bundles[-1].append(f)

fastpli.io.fiber.save('/tmp/fastpli.example.model_solver.dat', fiber_bundles)
''' setup solver
'''
solver = fastpli.model.Solver()
solver.fiber_bundles = filled_fiber_bundles
solver.drag = 0
solver.obj_min_radius = 0
solver.obj_mean_length = 15
solver.omp_num_threads = 8

solved_fbs = solver.fiber_bundles
fastpli.io.fiber.save('/tmp/fastpli.example.model_solver_.dat', solved_fbs)
''' run solver
'''
for i in range(1000):

    solved = solver.step()

    if i % 10 == 0:
        print("step:", i, solver.num_obj, solver.num_col_obj)
        solver.draw_scene()

    if solved:
        print("step:", i, solver.num_obj, solver.num_col_obj)
        solver.draw_scene()
        break
