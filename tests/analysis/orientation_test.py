import unittest

import numpy as np

import fastpli.io
import fastpli.analysis

np.random.seed(42)


class MainTest(unittest.TestCase):
    def test_remap_direction(self):
        phi = np.linspace(-42 * np.pi, 42 * np.pi, 1000)
        phi = fastpli.analysis.orientation.remap_direction(phi)
        self.assertTrue(np.all(phi >= 0) or np.all(phi < np.pi))

    def test_remap_orientation(self):
        phi, theta = np.mgrid[-42 * np.pi:42 * np.pi:100j,
                              -42 * np.pi:42 * np.pi:100j]
        phi = phi.ravel()
        theta = theta.ravel()
        phi, theta = fastpli.analysis.orientation.remap_orientation(phi, theta)
        x = np.multiply(np.cos(phi), np.sin(theta))
        y = np.multiply(np.sin(phi), np.sin(theta))
        z = np.cos(theta)
        x[z < 0] = -x[z < 0]
        y[z < 0] = -y[z < 0]
        z[z < 0] = -z[z < 0]
        phi_ = np.arctan2(y, x)
        theta_ = np.arccos(z)
        phi_[phi_ < 0] += 2 * np.pi
        self.assertTrue(np.allclose(phi, phi_))
        self.assertTrue(np.allclose(theta, theta_))

    def test_remap_spherical(self):
        phi, theta = np.mgrid[-42 * np.pi:42 * np.pi:100j,
                              -42 * np.pi:42 * np.pi:100j]
        phi = phi.ravel()
        theta = theta.ravel()
        phi, theta = fastpli.analysis.orientation.remap_spherical(phi, theta)
        x = np.multiply(np.cos(phi), np.sin(theta))
        y = np.multiply(np.sin(phi), np.sin(theta))
        z = np.cos(theta)
        phi_ = np.arctan2(y, x)
        theta_ = np.arccos(z)
        phi_[phi_ < 0] += 2 * np.pi
        self.assertTrue(np.allclose(phi, phi_))
        self.assertTrue(np.allclose(theta, theta_))

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
