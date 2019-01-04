import numpy as np
import fastpli

# setup fibers
VOLUME = 10
NFIBER = 100
np.random.seed(42)

# 2d plane with rnd seedpoints for fibers
planes = np.empty((NFIBER, 2, 3))
planes[:, :, 0:2] = np.random.uniform(-VOLUME, +VOLUME, (NFIBER, 2, 2))
planes[:, 0, 2] = -VOLUME
planes[:, 1, 2] = +VOLUME

fiber_bundles = [[]]
for points in planes:
    fiber_bundles[0].append(
        fastpli.objects.Fiber(points, np.random.uniform(0.5, 1.0, 2)))

# setup solver
solver = fastpli.model.Solver()
solver.set_fiber_bundles(fiber_bundles)
solver.set_parameters(drag=0, obj_min_radius=10, obj_mean_length=1)
# solver.set_col_voi([0, 0, 0], [10, 10, 10])
solver.set_omp_num_threads(8)

# run solver and plot results
# vis = fastpli.model.Visualizer()
for i in range(1000):
    print("step:", i, solver.num_obj, solver.num_col_obj)

    if solver.step():
        break

    if i % 25 == 0:
        # vis.set_fbs(solver.get_fiber_bundles())
        # vis.draw()
        solver.draw_scene()
