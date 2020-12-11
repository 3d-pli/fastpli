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

        phi, theta = fastpli.analysis.orientation.remap_orientation(
            np.deg2rad(0), np.deg2rad(0))
        self.assertTrue(phi == np.deg2rad(0))
        self.assertTrue(theta == np.deg2rad(0))

        phi, theta = fastpli.analysis.orientation.remap_orientation(
            np.deg2rad(90), np.deg2rad(10))
        self.assertTrue(phi == np.deg2rad(90))
        self.assertTrue(theta == np.deg2rad(10))

        phi, theta = fastpli.analysis.orientation.remap_orientation(
            np.deg2rad(360), np.deg2rad(25))
        self.assertTrue(phi == np.deg2rad(0))
        self.assertTrue(theta == np.deg2rad(25))

        phi, theta = fastpli.analysis.orientation.remap_orientation(
            np.deg2rad(25), np.deg2rad(123))
        self.assertAlmostEqual(phi, np.deg2rad(205))
        self.assertAlmostEqual(theta, (np.pi - np.deg2rad(123)))

        phi, theta = fastpli.analysis.orientation.remap_orientation(
            np.deg2rad(25), np.deg2rad(200))
        self.assertAlmostEqual(phi, np.deg2rad(205))
        self.assertAlmostEqual(theta, np.deg2rad(20))

        phi, theta = fastpli.analysis.orientation.remap_orientation(
            np.deg2rad(-25), np.deg2rad(-200))
        self.assertAlmostEqual(phi, np.deg2rad(180 - 25))
        self.assertAlmostEqual(theta, np.deg2rad(20))

        phi, theta = fastpli.analysis.orientation.remap_orientation(
            np.deg2rad(-25), np.deg2rad(200))
        self.assertAlmostEqual(phi, np.deg2rad(180 - 25))
        self.assertAlmostEqual(theta, np.deg2rad(20))

    def test_remap_spherical(self):
        phi = np.random.uniform(-1000, 1000, 1000)
        theta = np.random.uniform(-1000, 1000, 1000)
        phi, theta = fastpli.analysis.orientation.remap_spherical(phi, theta)
        self.assertTrue(
            np.all(phi >= 0) or np.all(phi < 2 * np.pi) or np.all(theta >= 0) or
            np.all(theta < np.pi))

        phi, theta = fastpli.analysis.orientation.remap_orientation(phi, theta)
        phi_, theta_ = fastpli.analysis.orientation.remap_spherical(phi, theta)
        self.assertTrue(
            np.array_equal(phi, phi_) and np.array_equal(theta, theta_))

        phi, theta = fastpli.analysis.orientation.remap_spherical(
            np.deg2rad(0), np.deg2rad(0))
        self.assertTrue(phi == np.deg2rad(0))
        self.assertTrue(theta == np.deg2rad(0))

        phi, theta = fastpli.analysis.orientation.remap_spherical(
            np.deg2rad(90), np.deg2rad(10))
        self.assertTrue(phi == np.deg2rad(90))
        self.assertTrue(theta == np.deg2rad(10))

        phi, theta = fastpli.analysis.orientation.remap_spherical(
            np.deg2rad(360), np.deg2rad(25))
        self.assertAlmostEqual(phi, np.deg2rad(0))
        self.assertAlmostEqual(theta, np.deg2rad(25))

        phi, theta = fastpli.analysis.orientation.remap_spherical(
            np.deg2rad(25), np.deg2rad(123))
        self.assertAlmostEqual(phi, np.deg2rad(25))
        self.assertAlmostEqual(theta, np.deg2rad(123))

        phi, theta = fastpli.analysis.orientation.remap_spherical(
            np.deg2rad(25), np.deg2rad(200))
        self.assertAlmostEqual(phi, np.deg2rad(180 + 25))
        self.assertAlmostEqual(theta, np.deg2rad(360 - 200))

        phi, theta = fastpli.analysis.orientation.remap_spherical(
            np.deg2rad(-25), np.deg2rad(-200))
        self.assertAlmostEqual(phi, np.deg2rad(360 - 25))
        self.assertAlmostEqual(theta, np.deg2rad(160))

        phi, theta = fastpli.analysis.orientation.remap_spherical(
            np.deg2rad(-25), np.deg2rad(200))
        self.assertAlmostEqual(phi, np.deg2rad((360 - 25) - 180))
        self.assertAlmostEqual(theta, np.deg2rad(360 - 200))

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
