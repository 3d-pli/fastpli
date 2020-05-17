import unittest

import matplotlib.pyplot as plt
import numpy as np

import fastpli.io
import fastpli.analysis


class MainTest(unittest.TestCase):

    def test_fiber_bundles(self):
        phi, theta = fastpli.analysis.orientation.fiber_bundles(
            fastpli.io.fiber_bundles.load("examples/cube.dat"))

        self.assertTrue(
            np.all(phi >= 0) or np.all(phi < 2 * np.pi) or np.all(theta >= 0) or
            np.all(theta < 0.5 * np.pi))

    def test_histogram(self):

        _, ax = plt.subplots(subplot_kw=dict(projection="polar"))
        phi = np.random.normal(np.pi / 3, 0.25, 100)
        theta = np.random.normal(np.pi / 4, 0.25, 100)
        _, _, _, pc = fastpli.analysis.orientation.histogram(phi, theta, ax)

        plt.colorbar(pc, ax=ax)
        ax.set_rmax(90)
        ax.set_rticks(range(0, 90, 10))
        ax.set_rlabel_position(22.5)
        ax.set_yticklabels([])
        ax.grid(True)
