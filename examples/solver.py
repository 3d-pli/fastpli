import fastpli.model.solver
import fastpli.model.sandbox
import fastpli.io

import numpy as np
import os

np.random.seed(42)

FILE_NAME = os.path.abspath(__file__)
FILE_PATH = os.path.dirname(FILE_NAME)
FILE_BASE = os.path.basename(FILE_NAME)
FILE_OUT = os.path.join(FILE_PATH, f'fastpli.example.{FILE_BASE}')

solver = fastpli.model.solver.Solver()
example = 'curved'

if example == 'curved':
    phi = np.linspace(0, 1 / 2 * np.pi, 16)
    fiber_bundle_trj = 200 * np.array([np.cos(phi), np.sin(phi), 0 * phi]).T
    population = fastpli.model.sandbox.seeds.triangular_circle(20, 5)
    fiber_radii = np.random.uniform(2.0, 10.0, population.shape[0])
    fiber_bundle = fastpli.model.sandbox.build.bundle(fiber_bundle_trj,
                                                      population, fiber_radii)
    solver.fiber_bundles = [fiber_bundle]

elif example == 'crossing':
    fiber_bundle_trj_0 = [[-150, 0, 0], [150, 0, 0]]
    fiber_bundle_trj_1 = [[0, -150, 0], [0, 150, 0]]

    population = fastpli.model.sandbox.seeds.triangular_circle(20, 5)

    fiber_radii = np.random.uniform(2.0, 10.0, population.shape[0])
    fiber_bundle_0 = fastpli.model.sandbox.build.bundle(fiber_bundle_trj_0,
                                                        population, fiber_radii)

    fiber_radii = np.random.uniform(2.0, 10.0, population.shape[0])
    fiber_bundle_1 = fastpli.model.sandbox.build.bundle(fiber_bundle_trj_1,
                                                        population, fiber_radii)

    solver.fiber_bundles = [fiber_bundle_0, fiber_bundle_1]
    solver.obj_min_radius = 10
    solver.obj_mean_length = 30

# additional parameter, 0 means disabled
# solver.drag = 0
# solver.obj_min_radius = 0
# solver.obj_mean_length = 0
# solver.omp_num_threads = 1

# run solver
solver.toggle_axis(True)
solver.draw_scene()
for i in range(1000):
    solved = solver.step()

    # calculate current overlap ratio
    overlap = solver.overlap / solver.num_col_obj if solver.num_col_obj else 0
    if i % 5 == 0:
        print(
            f'step: {i}, {solver.num_obj}/{solver.num_col_obj} {round(overlap * 100)}%'
        )
        solver.draw_scene()
        # solver.save_ppm(f'solver_{i:03}.ppm')  # save a ppm image

    if solved:
        print(f'solved: {i}, {solver.num_obj}/{solver.num_col_obj}')
        solver.draw_scene()
        break

fastpli.io.fiber_bundles.save(f'{FILE_OUT}.dat', solver.fiber_bundles)

print('Done')
