import fastpli.model.solver
import fastpli.model.sandbox
import fastpli.io

import numpy as np
import os, sys

FILE_NAME = os.path.abspath(__file__)
FILE_PATH = os.path.dirname(FILE_NAME)
FILE_BASE = os.path.basename(FILE_NAME)

np.random.seed(42)
### create fiber_bundle(s) ###
fiber_bundle_trj = []
for phi in range(0, 91, 5):
    fiber_bundle_trj.append(
        [200 * np.sin(np.deg2rad(phi)), 200 * np.cos(np.deg2rad(phi)), 0])
fiber_bundle_trj = np.array(fiber_bundle_trj)

population = fastpli.model.sandbox.seeds.triangular_circle(20, 5)

fiber_radii = np.random.uniform(2.0, 5.0, population.shape[0])
fiber_bundle = fastpli.model.sandbox.build.bundle(fiber_bundle_trj, population,
                                                  fiber_radii)

file_name = 'fastpli.example.' + FILE_BASE + '.dat'
print(f"creating file: {file_name}")

fastpli.io.fiber_bundles.save(file_name, [fiber_bundle])

### setup solver ###
solver = fastpli.model.solver.Solver()
solver.fiber_bundles = [fiber_bundle]
solver.drag = 0
solver.obj_min_radius = 0
solver.obj_mean_length = 15
solver.omp_num_threads = 1

### run solver ###
solver.draw_scene()
for i in range(1000):
    solved = solver.step()

    overlap = solver.overlap / solver.num_col_obj if solver.num_col_obj else 0
    if i % 10 == 0:
        print(
            f"step: {i}, {solver.num_obj}/{solver.num_col_obj} {round(overlap * 100)}%"
        )
        solver.draw_scene()
        # solver.save_ppm(f'solver_{i}.ppm')

    if solved:
        print(f"solved: {i}, {solver.num_obj}/{solver.num_col_obj}")
        solver.draw_scene()
        break

file_name = 'fastpli.example.' + FILE_BASE + '.solved.dat'
print(f"creating file: {file_name}")

fastpli.io.fiber_bundles.save(file_name, solver.fiber_bundles)

print("Done")
