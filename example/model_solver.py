import numpy as np
import fastpli

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# setup fibers
VOLUME = 10
NFIBER = 50
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


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# run solver and plot results
for i in range(100):
    solved = solver.step()
    if not solved:
        for fb in solver.get_fiber_bundles():
            for f in fb:
                p = f.points
                ax.plot(p[:, 0], p[:, 1], p[:, 2])
    else:
        break

    print("step:", i, solver.num_obj, solver.num_col_obj)
    ax.set_xlim(-1.0 * VOLUME, 1.0 * VOLUME)
    ax.set_ylim(-1.0 * VOLUME, 1.0 * VOLUME)
    ax.set_zlim(-1.0 * VOLUME, 1.0 * VOLUME)
    ax.set_aspect('equal')
    plt.pause(0.1)
    ax.cla()
