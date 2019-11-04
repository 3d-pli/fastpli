import fastpli.model.sandbox as sandbox

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

PLOT = False


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


def plot_fiber_bundle(fb):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for fiber in fb:
        ax.plot(fiber[:, 0], fiber[:, 1], fiber[:, 2])
    set_axes_equal(ax)
    plt.show()


# create fiber bundle along trajectory
seeds = sandbox.seeds.triangular_grid(a=42, b=42, spacing=4, center=True)
circ_seeds = sandbox.seeds.crop_circle(radius=21, seeds=seeds)
t = np.linspace(0, 2 * np.pi, 50, True)
traj = np.array((10 * t, np.cos(t), np.zeros(t.size))).T
fiber_bundle = sandbox.build.bundle(traj=traj,
                                    seeds=circ_seeds,
                                    radii=np.random.uniform(
                                        0.5, 0.8, circ_seeds.shape[0]),
                                    scale=2 + 0.5 * np.sin(t))
if PLOT:
    plot_fiber_bundle(fiber_bundle)

# create circular shaped triangular seeds
seeds = sandbox.seeds.triangular_grid(a=200, b=200, spacing=5, center=True)
fiber_bundle = sandbox.build.cylinder(
    p=(0, 80, 50),
    q=(40, 80, 100),
    r_in=20,
    r_out=40,
    seeds=seeds,
    radii=1,
    alpha=np.deg2rad(20),
    beta=np.deg2rad(160),
    mode='r'  # 'c', 'p'
)
if PLOT:
    plot_fiber_bundle(fiber_bundle)

# create circular shaped triangular seeds
p = np.array([0, 80, 50])
q = np.array([40, 180, 100])
d = np.max(np.abs(p - q)) * np.sqrt(3)
seeds = sandbox.seeds.triangular_grid(a=d, b=d, spacing=5, center=True)
fiber_bundle = sandbox.build.cuboid(p=p,
                                    q=q,
                                    phi=np.deg2rad(45),
                                    theta=np.deg2rad(90),
                                    seeds=seeds,
                                    radii=1)
if True:
    plot_fiber_bundle(fiber_bundle)
