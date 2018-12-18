import unittest
import numpy as np

import fastpli
import fastpli.tools


class MainTest(unittest.TestCase):

    def setUp(self):
        self.points = [0, 0, 0, 0, 0, 1, 0, 0, 2]
        self.radii = [1, 1, 1]
        self.fiber = fastpli.objects.Fiber(self.points, self.radii)

    def test_init(self):
        self.assertIsInstance(self.fiber.points, (np.ndarray))
        self.assertIsInstance(self.fiber.radii, (np.ndarray))
        self.assertTrue(np.array_equal(
            np.array(self.points, dtype=np.float32).reshape(-1, 3), self.fiber.points))
        self.assertTrue(
            np.array_equal(
                np.array(
                    self.radii,
                    dtype=np.float32),
                self.fiber.radii))

    # def test_cast_base_class(self):
    #     test = fastpli.objects._fiber_cpp._FiberCPP([],[])
    #     test.__class__ = fastpli.objects.Fiber
    #     print(repr(test))

    def test_rotation(self):
        points = self.fiber.points
        rot = fastpli.tools.rotations.rot_theta_phi(
            np.deg2rad(30), np.deg2rad(60))
        self.fiber.rotate_fiber(rot)
        self.assertFalse(np.array_equal(self.fiber.points, points))

        rot = fastpli.tools.rotations.rot_theta_phi(
            np.deg2rad(-30), np.deg2rad(60))
        self.fiber.rotate_fiber(rot)
        self.assertTrue(np.sum(self.fiber.points - points) < 1e-7)

    def test_rotate_around_point(self):
        self.fiber.rotate_around_point(fastpli.tools.rotations.rot_theta_phi(
            np.deg2rad(30), np.deg2rad(60)), [0, 0, 0])

    def test_translate(self):
        self.fiber.translate([0, 0, 0])

    def test_scale_points(self):
        points = self.fiber.points
        self.fiber.scale_points(1)
        self.assertTrue(np.array_equal(points, self.fiber.points))

        self.fiber.scale_points(2.5)
        self.assertTrue(np.array_equal(points * 2.5, self.fiber.points))

    def test_scale_radii(self):
        radii = self.fiber.radii
        self.fiber.scale_radii(1)
        self.assertTrue(np.array_equal(radii, self.fiber.radii))

        self.fiber.scale_radii(2.5)
        self.assertTrue(np.array_equal(radii * 2.5, self.fiber.radii))

    def test_scale(self):
        (points, radii) = self.fiber.data
        self.fiber.scale(1)
        self.assertTrue(np.array_equal(points, self.fiber.points))
        self.assertTrue(np.array_equal(radii, self.fiber.radii))

        self.fiber.scale(2.5)
        self.assertTrue(np.array_equal(points * 2.5, self.fiber.points))
        self.assertTrue(np.array_equal(radii * 2.5, self.fiber.radii))


if __name__ == '__main__':
    unittest.main()
