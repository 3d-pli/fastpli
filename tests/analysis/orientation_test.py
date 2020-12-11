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
        phi = np.linspace(-42 * np.pi, 42 * np.pi, 1000)
        theta = np.linspace(-42 * np.pi, 42 * np.pi, 1000)
        phi, theta = fastpli.analysis.orientation.remap_orientation(phi, theta)
        self.assertTrue(
            np.all(phi >= 0) or np.all(phi < 2 * np.pi) or np.all(theta >= 0) or
            np.all(theta <= 0.5 * np.pi))

        phi, theta = fastpli.analysis.orientation.remap_orientation(0, 0)
        self.assertTrue(phi == 0 and theta == 0, f"{phi}, {theta}")

        phi, theta = fastpli.analysis.orientation.remap_orientation(
            0, 0.5 * np.pi)
        self.assertTrue(phi == 0 and theta == 0.5 * np.pi, f"{phi}, {theta}")

        phi, theta = fastpli.analysis.orientation.remap_orientation(
            2 * np.pi, 0)
        self.assertTrue(phi == 0 and theta == 0, f"{phi}, {theta}")

        phi, theta = fastpli.analysis.orientation.remap_orientation(
            2 * np.pi, 0.5 * np.pi)
        self.assertTrue(phi == 0 and theta == 0.5 * np.pi, f"{phi}, {theta}")

        phi, theta = fastpli.analysis.orientation.remap_orientation(
            0, 0.75 * np.pi)
        self.assertTrue(phi == np.pi and theta == 0.25 * np.pi,
                        f"{phi}, {theta}")

        phi, theta = fastpli.analysis.orientation.remap_orientation(
            0, -0.75 * np.pi)
        self.assertTrue(phi == 0 and theta == 0.25 * np.pi, f"{phi}, {theta}")

        phi, theta = fastpli.analysis.orientation.remap_orientation(
            -0.5 * np.pi, 0)
        self.assertTrue(phi == 1.5 * np.pi and theta == 0, f"{phi}, {theta}")

    def test_remap_spherical(self):
        phi = np.linspace(-42 * np.pi, 42 * np.pi, 1000)
        theta = np.linspace(-42 * np.pi, 42 * np.pi, 1000)
        phi, theta = fastpli.analysis.orientation.remap_spherical(phi, theta)
        self.assertTrue(
            np.all(phi >= 0) or np.all(phi < 2 * np.pi) or np.all(theta >= 0) or
            np.all(theta < np.pi))

        phi, theta = fastpli.analysis.orientation.remap_orientation(phi, theta)
        phi_, theta_ = fastpli.analysis.orientation.remap_spherical(phi, theta)
        self.assertTrue(
            np.array_equal(phi, phi_) and np.array_equal(theta, theta_))

        phi, theta = fastpli.analysis.orientation.remap_spherical(0, 0)
        self.assertTrue(phi == 0 and theta == 0, f"{phi}, {theta}")

        phi, theta = fastpli.analysis.orientation.remap_spherical(
            0, 0.5 * np.pi)
        self.assertTrue(phi == 0 and theta == 0.5 * np.pi, f"{phi}, {theta}")

        phi, theta = fastpli.analysis.orientation.remap_spherical(0, np.pi)
        self.assertTrue(phi == 0 and theta == np.pi, f"{phi}, {theta}")

        phi, theta = fastpli.analysis.orientation.remap_spherical(2 * np.pi, 0)
        self.assertTrue(phi == 0 and theta == 0, f"{phi}, {theta}")

        phi, theta = fastpli.analysis.orientation.remap_spherical(
            2 * np.pi, 0.5 * np.pi)
        self.assertTrue(phi == 0 and theta == 0.5 * np.pi, f"{phi}, {theta}")

        phi, theta = fastpli.analysis.orientation.remap_spherical(
            0, 0.75 * np.pi)
        self.assertTrue(phi == 0 and theta == 0.75 * np.pi, f"{phi}, {theta}")

        phi, theta = fastpli.analysis.orientation.remap_spherical(
            0, -0.75 * np.pi)
        self.assertTrue(phi == np.pi and theta == 0.75 * np.pi,
                        f"{phi}, {theta}")

        phi, theta = fastpli.analysis.orientation.remap_spherical(
            3.25 * np.pi, 1.75 * np.pi)
        self.assertTrue(
            np.isclose(phi, 0.25 * np.pi) and np.isclose(theta, 0.25 * np.pi),
            f"{phi}, {theta}")

        phi, theta = fastpli.analysis.orientation.remap_spherical(
            3.25 * np.pi, -1.75 * np.pi)
        self.assertTrue(
            np.isclose(phi, 1.25 * np.pi) and np.isclose(theta, 0.25 * np.pi),
            f"{phi}, {theta}")

        phi, theta = fastpli.analysis.orientation.remap_spherical(
            -0.5 * np.pi, 0)
        self.assertTrue(phi == 1.5 * np.pi and theta == 0, f"{phi}, {theta}")

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
