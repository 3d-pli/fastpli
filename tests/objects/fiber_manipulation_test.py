import unittest
import numpy as np

import fastpli.objects
import fastpli.tools


class MainTest(unittest.TestCase):

    # TODO: implement object.fiber.*manipulations*
    def setUp(self):
        self.fiber = np.array([[0, 0, 0, 1], [1, 1, 1, 2]], dtype=float)
        self.fiber_bundle = [self.fiber.copy()]
        self.fiber_bundles = [[self.fiber.copy()]]

    def test_resize(self):
        fiber = fastpli.objects.fiber.Rescale(self.fiber, 10)
        self.assertTrue(np.array_equal(fiber, self.fiber * 10))

        fb = fastpli.objects.fiber_bundle.Rescale(self.fiber_bundle, 10)
        for f in fb:
            self.assertTrue(np.array_equal(f, self.fiber * 10))

        fbs = fastpli.objects.fiber_bundles.Rescale(self.fiber_bundles, 10)
        for fb in fbs:
            for f in fb:
                self.assertTrue(np.array_equal(f, self.fiber * 10))

        fiber = fastpli.objects.fiber.Rescale(self.fiber, 10, mod='points')
        self.assertTrue(np.array_equal(fiber[:, :-2], self.fiber[:, :-2] * 10))
        self.assertTrue(np.array_equal(fiber[:, -1], self.fiber[:, -1]))

        fiber = fastpli.objects.fiber.Rescale(self.fiber, 10, mod='radii')
        self.assertTrue(np.array_equal(fiber[:, :-2], self.fiber[:, :-2]))
        self.assertTrue(np.array_equal(fiber[:, -1], self.fiber[:, -1] * 10))

    def test_rotation(self):
        fiber = fastpli.objects.fiber.Rotate(self.fiber,
                                             fastpli.tools.rotation.x(0))
        self.assertTrue(np.array_equal(self.fiber, fiber))

        fiber = fastpli.objects.fiber.Rotate(
            self.fiber, fastpli.tools.rotation.x(np.deg2rad(90)))
        self.assertTrue(
            np.allclose(fiber, np.array([[0, 0, 0, 1], [1, -1, 1, 2]])))

        fiber = fastpli.objects.fiber.Rotate(
            self.fiber, fastpli.tools.rotation.x(np.deg2rad(90)), [1, 1, 1])
        self.assertTrue(
            np.allclose(fiber, np.array([[0, 2, 0, 1], [1, 1, 1, 2]])))

        for f in self.fiber_bundle:
            fiber = fastpli.objects.fiber.Rotate(
                f, fastpli.tools.rotation.x(np.deg2rad(90)), [1, 1, 1])
            self.assertTrue(
                np.allclose(fiber, np.array([[0, 2, 0, 1], [1, 1, 1, 2]])))

        for fb in self.fiber_bundles:
            for f in fb:
                fiber = fastpli.objects.fiber.Rotate(
                    f, fastpli.tools.rotation.x(np.deg2rad(90)), [1, 1, 1])
                self.assertTrue(
                    np.allclose(fiber, np.array([[0, 2, 0, 1], [1, 1, 1, 2]])))

    def test_translate(self):
        fiber = fastpli.objects.fiber.Translate(self.fiber, [1, 1, 1])
        self.assertTrue(
            np.array_equal(fiber[:, :3],
                           self.fiber[:, :3] + np.array([1, 1, 1])))
        self.assertTrue(np.array_equal(fiber[:, -1], self.fiber[:, -1]))

        for f in self.fiber_bundle:
            fiber = fastpli.objects.fiber.Translate(f, [1, 1, 1])
            self.assertTrue(
                np.array_equal(fiber[:, :3],
                               self.fiber[:, :3] + np.array([1, 1, 1])))
            self.assertTrue(np.array_equal(f[:, -1], self.fiber[:, -1]))

        for fb in self.fiber_bundles:
            for f in fb:
                fiber = fastpli.objects.fiber.Translate(f, [1, 1, 1])
                self.assertTrue(
                    np.array_equal(fiber[:, :3],
                                   self.fiber[:, :3] + np.array([1, 1, 1])))
                self.assertTrue(np.array_equal(f[:, -1], self.fiber[:, -1]))


if __name__ == '__main__':
    unittest.main()
