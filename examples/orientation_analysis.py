import fastpli.analysis
import fastpli.io
import fastpli.tools

import numpy as np
import os
import matplotlib.pyplot as plt

np.random.seed(42)

FILE_NAME = os.path.abspath(__file__)
FILE_PATH = os.path.dirname(FILE_NAME)
FILE_BASE = os.path.basename(FILE_NAME)

# fiber_bundles = fastpli.io.fiber_bundles.load(
#     os.path.join(FILE_PATH, 'cube.dat'))

fiber_bundles = [[]]

OMEGA = 20
THETA = 42
PHI = 42
N = 4200

# rnd orientation along z
theta = np.random.normal(0, np.deg2rad(OMEGA), N)
phi = np.random.normal(0, np.deg2rad(360), N)

for p, t in zip(phi, theta):
    v = np.array([0, 0, 1])
    rot = fastpli.tools.rotation.y(t)
    v = np.dot(rot, v)
    rot = fastpli.tools.rotation.z(p)
    v = np.dot(rot, v)

    # rotate into main orientation
    rot = fastpli.tools.rotation.y(np.deg2rad(THETA))
    v = np.dot(rot, v)
    rot = fastpli.tools.rotation.z(np.deg2rad(PHI))
    v = np.dot(rot, v)

    fiber_bundles[0].append(np.vstack([[0, 0, 0], v]))

phi, theta = fastpli.analysis.orientation.fiber_bundles(fiber_bundles)

_, ax = plt.subplots(subplot_kw=dict(projection='polar'))
_, _, _, pc = fastpli.analysis.orientation.histogram(phi,
                                                     theta,
                                                     ax=ax,
                                                     n_phi=60,
                                                     n_theta=30,
                                                     weight_area=False)
cbar = plt.colorbar(pc, ax=ax)
cbar.ax.set_title('#')
ax.set_rmax(90)
ax.set_rticks(range(0, 90, 10))
ax.set_rlabel_position(22.5)
ax.set_yticklabels([])
ax.grid(True)

plt.show()
