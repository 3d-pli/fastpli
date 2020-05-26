import unittest

import numpy as np

import fastpli.io
import fastpli.analysis

np.random.seed(42)


class MainTest(unittest.TestCase):

    def test_remap_direction(self):
        phi = np.random.uniform(-1000, 1000, 1000)
        phi = fastpli.analysis.orientation.remap_direction(phi)
        self.assertTrue(np.all(phi >= 0) or np.all(phi < np.pi))

    def test_remap_orientation(self):
        phi = np.random.uniform(-1000, 1000, 1000)
        theta = np.random.uniform(-1000, 1000, 1000)
        phi, theta = fastpli.analysis.orientation.remap_orientation(phi, theta)
        self.assertTrue(
            np.all(phi >= 0) or np.all(phi < 2 * np.pi) or np.all(theta >= 0) or
            np.all(theta < 0.5 * np.pi))

    def test_fiber_bundles(self):
        fastpli.analysis.orientation.fiber_bundles(
            fastpli.io.fiber_bundles.load("examples/cube.dat"))

    def test_histogram(self):
        phi = np.random.normal(np.pi / 3, 0.5, 1000)
        theta = np.random.normal(np.deg2rad(45), 0.5, 1000)

        fastpli.analysis.orientation.histogram(phi,
                                               theta,
                                               n_phi=60,
                                               n_theta=30,
                                               weight_area=True)

        fastpli.analysis.orientation.histogram(phi,
                                               theta,
                                               n_phi=60,
                                               n_theta=30,
                                               weight_area=False)
