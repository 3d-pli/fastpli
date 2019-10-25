import fastpli.model.sandbox

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def set_axes_equal(ax):
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


# set parameter
p0 = (0, 80, 50)
p1 = (40, 80, 100)
r0 = 20
r1 = 40
alpha = np.deg2rad(20)
beta = alpha + np.deg2rad(180)
mode = 'r'  # 'p' 'c' 'r'
spacing = 5
radius = 1

# create circular shaped triangular seeds
seeds = fastpli.model.sandbox.seeds.triangle_grid(400, 400, spacing)
seeds = seeds - np.array([200, 200])

fiber_bundle = fastpli.model.sandbox.build.cylinder(p0, p1, r0, r1, seeds,
                                                    radius, alpha, beta, mode)

# plot trajectories
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for fiber in fiber_bundle:
    ax.plot(fiber[:, 0], fiber[:, 1], fiber[:, 2])
set_axes_equal(ax)
plt.show()
