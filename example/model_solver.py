import numpy as np
import fastpli

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


fiber_bundles = [[
    fastpli.objects.Fiber([-1, -1, -1, 1, 1, 1], [1, 1]),
    fastpli.objects.Fiber([1, 1, -1, -1, -1, 1], [1, 1])]]

solver = fastpli.model.Solver()

solver.set_fibers(fiber_bundles)
solver.set_parameter(drag=0, obj_min_radius=10, obj_mean_length=1)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for i in range(100):
    solved = solver.step()
    print(i, solved)

    if not solved:
        print("num:", solver.num_objects)
        for fb in solver.get_fibers():
            for f in fb:
                p = f.points
                ax.plot(p[:, 0], p[:, 1], p[:, 2])
    else:
        break

    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    ax.set_zlim(-2, 2)
    plt.pause(0.1)
    ax.cla()
    print()
