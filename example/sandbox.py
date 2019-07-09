import fastpli.model.sandbox
import fastpli.io
import fastpli.tools

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# ax.set_aspect('equal')


def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


p0 = (0, 0, 0)
p1 = (100, 100, 100)
r0 = 20
r1 = 40
alpha = np.deg2rad(22.5)
beta = alpha + np.deg2rad(90)
mode = 'radial'
spacing = 5
steps = 100
radius = 1

data = fastpli.model.sandbox.shape.cylinder(p0, p1, r0, r1, alpha, beta, mode,
                                            spacing, radius, steps)

for d in data:
    ax.plot(d[:, 0], d[:, 1], d[:, 2])

set_axes_equal(ax)
plt.show()
