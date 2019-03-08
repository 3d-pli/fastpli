import fastpli
from fastpli.model import bundle_generator

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
        [200 * np.sin(np.deg2rad(phi)),
         200 * np.cos(np.deg2rad(phi)),
         0,
         1])
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

''' setup solver
'''
solver = fastpli.model.Solver()
solver.fiber_bundles = filled_fiber_bundles
solver.drag = 0
solver.obj_min_radius = 0
solver.obj_mean_length = 15
solver.omp_num_threads = 8


# solved_fbs = solver.fiber_bundles
# with open(os.path.join(FILE_PATH,'test.dat'), 'w') as file:
#     for fb in solved_fbs:
#         for f in fb:
#             p,r = f.data
#             for i in range(r.size):
# file.write(str(p[i,0]) + '\t' + str(p[i,1]) + '\t' + str(p[i,2]) + '\t'
# + str(r[i]) + '\n')

#             file.write('\n')


''' run solver
'''
for i in range(1000):
    solver.step()

    if i % 10 == 0:
        print("step:", i, solver.num_obj, solver.num_col_obj)
        solver.draw_scene()

# solver.draw_scene()
# while(True):
#     time.sleep(0.1)


# solved_fbs = solver.fiber_bundles
# with open(os.path.join(FILE_PATH,'test.dat'), 'w') as file:
#     for fb in solved_fbs:
#         for f in fb:
#             p,r = f.data
#             for i in range(r.size):
# file.write(str(p[i,0]) + '\t' + str(p[i,1]) + '\t' + str(p[i,2]) + '\t'
# + str(r[i]) + '\n')

#             file.write('\n')
