import unittest

import numpy as np

import fastpli.io
import fastpli.analysis

np.random.seed(42)


class MainTest(unittest.TestCase):

    def test_remap_direction(self):
        phi = fastpli.analysis.orientation.remap_direction(0)
        self.assertTrue(phi == 0)

        phi = np.linspace(-42 * np.pi, 42 * np.pi, 1000)
        phi = fastpli.analysis.orientation.remap_direction(phi)
        self.assertTrue(np.all(phi >= 0) or np.all(phi < np.pi))

    def test_remap_sphere(self):

        phi, theta = fastpli.analysis.orientation.remap_sphere(0, 0)
        self.assertTrue(phi == 0)
        self.assertTrue(theta == 0)
        phi, theta = fastpli.analysis.orientation.remap_sphere(0, np.pi)
        self.assertTrue(phi == 0)
        self.assertTrue(theta == np.pi)
        phi, theta = fastpli.analysis.orientation.remap_sphere(2 * np.pi, 0)
        self.assertTrue(phi == 0)
        self.assertTrue(theta == 0)
        phi, theta = fastpli.analysis.orientation.remap_sphere(2 * np.pi, np.pi)
        self.assertTrue(phi == 0)
        self.assertTrue(theta == np.pi)

        # float64
        phi, theta = np.mgrid[-4 * np.pi:4 * np.pi:121j,
                              -4 * np.pi:4 * np.pi:121j]
        phi = phi.ravel()
        theta = theta.ravel()
        phi, theta = fastpli.analysis.orientation.remap_sphere(phi, theta)
        self.assertTrue(np.all(phi >= 0))
        self.assertTrue(np.all(phi < 2 * np.pi))
        self.assertTrue(np.all(theta >= 0))
        self.assertTrue(np.all(theta <= np.pi))

        x = np.multiply(np.cos(phi), np.sin(theta))
        y = np.multiply(np.sin(phi), np.sin(theta))
        z = np.cos(theta)
        phi_ = np.arctan2(y, x)
        theta_ = np.arccos(z)
        phi_[phi_ < 0] += 2 * np.pi
        self.assertTrue(np.allclose(phi, phi_))
        self.assertTrue(np.allclose(theta, theta_))

        # float32
        phi, theta = np.mgrid[-4 * np.pi:4 * np.pi:121j,
                              -4 * np.pi:4 * np.pi:121j]
        phi = phi.ravel().astype(np.float32)
        theta = theta.ravel().astype(np.float32)
        phi, theta = fastpli.analysis.orientation.remap_sphere(phi, theta)
        x = np.multiply(np.cos(phi), np.sin(theta))
        y = np.multiply(np.sin(phi), np.sin(theta))
        z = np.cos(theta)
        phi_ = np.arctan2(y, x)
        theta_ = np.arccos(z)
        phi_[phi_ < 0] += 2 * np.pi
        self.assertTrue(np.allclose(phi, phi_, atol=1e-6))
        self.assertTrue(np.allclose(theta, theta_, atol=1e-6))

    def test_remap_half_sphere_z(self):
        phi, theta = fastpli.analysis.orientation.remap_half_sphere_z(0, 0)
        self.assertTrue(phi == 0)
        self.assertTrue(theta == 0)
        phi, theta = fastpli.analysis.orientation.remap_half_sphere_z(
            2 * np.pi, 0)
        self.assertTrue(phi == 0)
        self.assertTrue(theta == 0)
        phi, theta = fastpli.analysis.orientation.remap_half_sphere_z(
            0, np.pi / 2)
        self.assertTrue(phi == 0)
        self.assertTrue(theta == np.pi / 2)
        phi, theta = fastpli.analysis.orientation.remap_half_sphere_z(
            2 * np.pi, np.pi / 2)
        self.assertTrue(phi == 0)
        self.assertTrue(theta == np.pi / 2)
        phi, theta = fastpli.analysis.orientation.remap_half_sphere_z(
            3 * np.pi, 3 * np.pi / 2)
        self.assertTrue(phi == 0)
        self.assertTrue(theta == np.pi / 2)

        # float64
        phi, theta = np.mgrid[-4 * np.pi:4 * np.pi:121j,
                              -4 * np.pi:4 * np.pi:121j]
        phi = phi.ravel()
        theta = theta.ravel()
        phi, theta = fastpli.analysis.orientation.remap_half_sphere_z(
            phi, theta)

        self.assertTrue(np.all(phi >= 0))
        self.assertTrue(np.all(phi < 2 * np.pi))
        self.assertTrue(np.all(theta >= 0))
        self.assertTrue(np.all(theta < np.pi))

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

        # float32
        phi, theta = np.mgrid[-4 * np.pi:4 * np.pi:121j,
                              -4 * np.pi:4 * np.pi:121j]
        phi = phi.ravel().astype(np.float32)
        theta = theta.ravel().astype(np.float32)
        phi, theta = fastpli.analysis.orientation.remap_half_sphere_z(
            phi, theta)
        x = np.multiply(np.cos(phi), np.sin(theta))
        y = np.multiply(np.sin(phi), np.sin(theta))
        z = np.cos(theta)
        x[z < 0] = -x[z < 0]
        y[z < 0] = -y[z < 0]
        z[z < 0] = -z[z < 0]
        phi_ = np.arctan2(y, x)
        theta_ = np.arccos(z)
        phi_[phi_ < 0] += 2 * np.pi
        self.assertTrue(np.allclose(phi, phi_, atol=1e-6))
        self.assertTrue(np.allclose(theta, theta_, atol=1e-6))

    def test_remap_half_sphere_x(self):

        phi, theta = fastpli.analysis.orientation.remap_half_sphere_x(0, 0)
        self.assertTrue(phi == 0)
        self.assertTrue(theta == 0)
        phi, theta = fastpli.analysis.orientation.remap_half_sphere_x(
            2 * np.pi, 0)
        self.assertTrue(phi == 0)
        self.assertTrue(theta == 0)
        phi, theta = fastpli.analysis.orientation.remap_half_sphere_x(
            0, np.pi / 2)
        self.assertTrue(phi == 0)
        self.assertTrue(theta == np.pi / 2)
        phi, theta = fastpli.analysis.orientation.remap_half_sphere_x(
            2 * np.pi, np.pi / 2)
        self.assertTrue(phi == 0)
        self.assertTrue(theta == np.pi / 2)
        phi, theta = fastpli.analysis.orientation.remap_half_sphere_x(
            3 * np.pi, 3 * np.pi / 2)
        self.assertTrue(phi == 0)
        self.assertTrue(theta == np.pi / 2)
        phi, theta = fastpli.analysis.orientation.remap_half_sphere_x(
            3 * np.pi / 2, 3 * np.pi / 2)
        self.assertTrue(phi == -np.pi / 2)
        self.assertTrue(theta == np.pi / 2)

        # float64
        phi, theta = np.mgrid[-4 * np.pi:4 * np.pi:121j,
                              -4 * np.pi:4 * np.pi:121j]
        phi = phi.ravel()
        theta = theta.ravel()
        phi, theta = fastpli.analysis.orientation.remap_half_sphere_x(
            phi, theta)

        self.assertTrue(np.all(phi >= -np.pi / 2))
        self.assertTrue(np.all(phi < np.pi / 2))
        self.assertTrue(np.all(theta >= 0))
        self.assertTrue(np.all(theta < np.pi))

        x = np.multiply(np.cos(phi), np.sin(theta))
        y = np.multiply(np.sin(phi), np.sin(theta))
        z = np.cos(theta)
        y[x < 0] = -y[x < 0]
        z[x < 0] = -z[x < 0]
        x[x < 0] = -x[x < 0]
        phi_ = np.arctan2(y, x)
        theta_ = np.arccos(z)
        self.assertTrue(np.allclose(phi, phi_))
        self.assertTrue(np.allclose(theta, theta_))

        # float64
        phi, theta = np.mgrid[-4 * np.pi:4 * np.pi:121j,
                              -4 * np.pi:4 * np.pi:121j]
        phi = phi.ravel().astype(np.float32)
        theta = theta.ravel().astype(np.float32)
        phi, theta = fastpli.analysis.orientation.remap_half_sphere_x(
            phi, theta)

        x = np.multiply(np.cos(phi), np.sin(theta))
        y = np.multiply(np.sin(phi), np.sin(theta))
        z = np.cos(theta)
        y[x < 0] = -y[x < 0]
        z[x < 0] = -z[x < 0]
        x[x < 0] = -x[x < 0]
        phi_ = np.arctan2(y, x)
        theta_ = np.arccos(z)
        self.assertTrue(np.allclose(phi, phi_, atol=1e-6))
        self.assertTrue(np.allclose(theta, theta_, atol=1e-6))

    def test_fiber_bundles(self):
        fastpli.analysis.orientation.fiber_bundles(
            fastpli.io.fiber_bundles.load('tests/cube.dat'))

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
