import unittest

import numpy as np

import fastpli.io
import fastpli.analysis

np.random.seed(42)


class MainTest(unittest.TestCase):

    def setUp(self):
        phi, theta = np.mgrid[-4 * np.pi:4 * np.pi:121j,
                              -4 * np.pi:4 * np.pi:121j]
        self.phi = phi.ravel()
        self.theta = theta.ravel()

    def test_remap_direction(self):
        phi = fastpli.analysis.orientation.remap_direction(0)
        self.assertTrue(phi == 0, f'{phi:.2f}')

        phi = np.linspace(-4 * np.pi, 4 * np.pi, 101)
        phi = fastpli.analysis.orientation.remap_direction(phi)
        self.assertTrue(np.all(phi >= 0) or np.all(phi < np.pi))

    def test_remap_sphere(self):
        phi, theta = fastpli.analysis.orientation.remap_sphere(0, 0)
        self.assertTrue(phi == 0 and theta == 0, f'{phi:.2f}, {theta:.2f}')
        phi, theta = fastpli.analysis.orientation.remap_sphere(0, np.pi)
        self.assertTrue(phi == 0 and theta == np.pi, f'{phi:.2f}, {theta:.2f}')
        phi, theta = fastpli.analysis.orientation.remap_sphere(2 * np.pi, 0)
        self.assertTrue(phi == 0 and theta == 0, f'{phi:.2f}, {theta:.2f}')
        phi, theta = fastpli.analysis.orientation.remap_sphere(2 * np.pi, np.pi)
        self.assertTrue(phi == 0 and theta == np.pi, f'{phi:.2f}, {theta:.2f}')

        # grid test
        phi, theta = fastpli.analysis.orientation.remap_sphere(
            self.phi, self.theta)
        self.assertTrue(np.all(phi >= 0))
        self.assertTrue(np.all(phi < 2 * np.pi))
        self.assertTrue(np.all(theta >= 0))
        self.assertTrue(np.all(theta <= np.pi))

        x = np.multiply(np.cos(self.phi), np.sin(self.theta))
        y = np.multiply(np.sin(self.phi), np.sin(self.theta))
        z = np.cos(self.theta)
        phi_ = np.arctan2(y, x)
        theta_ = np.arccos(z)
        phi_[phi_ < 0] += 2 * np.pi

        # edge cases
        phi_[theta_ == 0] = 0
        phi_[theta_ == np.pi] = 0
        phi_[np.isclose(phi_, 2 * np.pi)] = 0

        self.assertTrue(np.allclose(phi, phi_))
        self.assertTrue(np.allclose(theta, theta_))

    def test_remap_half_sphere_z(self):
        phi, theta = fastpli.analysis.orientation.remap_half_sphere_z(0, 0)
        self.assertTrue(phi == 0 and theta == 0, f'{phi:.2f}, {theta:.2f}')
        phi, theta = fastpli.analysis.orientation.remap_half_sphere_z(
            2 * np.pi, 0)
        self.assertTrue(phi == 0 and theta == 0, f'{phi:.2f}, {theta:.2f}')
        phi, theta = fastpli.analysis.orientation.remap_half_sphere_z(
            0, np.pi / 2)
        self.assertTrue(phi == 0 and theta == np.pi / 2,
                        f'{phi:.2f}, {theta:.2f}')
        phi, theta = fastpli.analysis.orientation.remap_half_sphere_z(
            2 * np.pi, np.pi / 2)
        self.assertTrue(phi == 0 and theta == np.pi / 2,
                        f'{phi:.2f}, {theta:.2f}')
        phi, theta = fastpli.analysis.orientation.remap_half_sphere_z(
            3 * np.pi, 3 * np.pi / 2)
        self.assertTrue(phi == 0 and theta == np.pi / 2,
                        f'{phi:.2f}, {theta:.2f}')
        phi, theta = fastpli.analysis.orientation.remap_half_sphere_z(
            0, 5 / 4 * np.pi / 2)
        self.assertTrue(phi == np.pi and theta == 3 / 4 * np.pi / 2,
                        f'{phi:.2f}, {theta:.2f}')

        # grid test
        phi, theta = fastpli.analysis.orientation.remap_half_sphere_z(
            self.phi, self.theta)

        self.assertTrue(np.all(phi >= 0))
        self.assertTrue(np.all(phi < 2 * np.pi))
        self.assertTrue(np.all(theta >= 0))
        self.assertTrue(np.all(theta < np.pi))

        x = np.multiply(np.cos(self.phi), np.sin(self.theta))
        y = np.multiply(np.sin(self.phi), np.sin(self.theta))
        z = np.cos(self.theta)
        x[z < 0] = -x[z < 0]
        y[z < 0] = -y[z < 0]
        z[z < 0] = -z[z < 0]
        phi_ = np.arctan2(y, x)
        theta_ = np.arccos(z)
        phi_[phi_ < 0] += 2 * np.pi

        # edge cases
        phi_[theta_ == 0] = 0
        phi_[np.logical_and(theta_ == np.pi / 2, phi_ >= np.pi)] -= np.pi
        phi_[np.isclose(phi_, 2 * np.pi)] = 0

        self.assertTrue(np.allclose(phi, phi_))
        self.assertTrue(np.allclose(theta, theta_))

    def test_remap_half_sphere_x(self):

        phi, theta = fastpli.analysis.orientation.remap_half_sphere_x(0, 0)
        self.assertTrue(phi == 0 and theta == 0, f'{phi:.2f}, {theta:.2f}')
        phi, theta = fastpli.analysis.orientation.remap_half_sphere_x(
            2 * np.pi, 0)
        self.assertTrue(phi == 0 and theta == 0, f'{phi:.2f}, {theta:.2f}')
        phi, theta = fastpli.analysis.orientation.remap_half_sphere_x(
            0, np.pi / 2)
        self.assertTrue(phi == 0 and theta == np.pi / 2,
                        f'{phi:.2f}, {theta:.2f}')
        phi, theta = fastpli.analysis.orientation.remap_half_sphere_x(
            2 * np.pi, np.pi / 2)
        self.assertTrue(phi == 0 and theta == np.pi / 2,
                        f'{phi:.2f}, {theta:.2f}')
        phi, theta = fastpli.analysis.orientation.remap_half_sphere_x(
            3 * np.pi, 3 * np.pi / 2)
        self.assertTrue(phi == 0 and theta == np.pi / 2,
                        f'{phi:.2f}, {theta:.2f}')
        phi, theta = fastpli.analysis.orientation.remap_half_sphere_x(
            3 * np.pi / 2, 3 * np.pi / 2)
        self.assertTrue(phi == -np.pi / 2 and theta == np.pi / 2,
                        f'{phi:.2f}, {theta:.2f}')
        phi, theta = fastpli.analysis.orientation.remap_half_sphere_x(
            0, 5 / 4 * np.pi / 2)
        self.assertTrue(phi == 0 and theta == 5 / 4 * np.pi / 2,
                        f'{phi:.2f}, {theta:.2f}')

        # grid test
        phi, theta = fastpli.analysis.orientation.remap_half_sphere_x(
            self.phi, self.theta)

        self.assertTrue(np.all(phi >= -np.pi / 2))
        self.assertTrue(np.all(phi < np.pi / 2))
        self.assertTrue(np.all(theta >= 0))
        self.assertTrue(np.all(theta < np.pi))

        x = np.multiply(np.cos(self.phi), np.sin(self.theta))
        y = np.multiply(np.sin(self.phi), np.sin(self.theta))
        z = np.cos(self.theta)
        y[x < 0] = -y[x < 0]
        z[x < 0] = -z[x < 0]
        x[x < 0] = -x[x < 0]
        phi_ = np.arctan2(y, x)
        theta_ = np.arccos(z)

        # edge cases
        phi_[theta_ == 0] = 0
        phi_[theta_ == np.pi] = 0
        theta_[np.isclose(theta_, np.pi)] = 0

        self.assertTrue(np.allclose(phi, phi_))
        self.assertTrue(np.allclose(theta, theta_))

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
