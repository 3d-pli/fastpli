import unittest
import numpy as np

import fastpli.objects
import fastpli.tools


class MainTest(unittest.TestCase):

    def setUp(self):
        self.fiber = fastpli.objects.Fiber([[0, 0, 0, 1], [1, 1, 1, 2]])
        self.fiber_bundle = fastpli.objects.FiberBundle(self.fiber.copy())
        self.fiber_bundles = fastpli.objects.FiberBundles(self.fiber.copy())

    def test_init(self):
        fastpli.objects.FiberBundle()
        fastpli.objects.FiberBundles()

        a = np.array([0, 0, 0, 0])
        _ = fastpli.objects.Fiber([[0, 0, 0, 1], [0, 0, 1, 2]])

        f = fastpli.objects.Fiber(a)
        self.assertTrue(isinstance(f, fastpli.objects.Fiber))
        f = fastpli.objects.Fiber(f)
        self.assertTrue(isinstance(f, fastpli.objects.Fiber))

        fb = fastpli.objects.FiberBundle([a])
        self.assertTrue(isinstance(fb, fastpli.objects.FiberBundle))
        fb = fastpli.objects.FiberBundle(f)
        self.assertTrue(isinstance(fb, fastpli.objects.FiberBundle))
        fb = fastpli.objects.FiberBundle(fb)
        self.assertTrue(isinstance(fb, fastpli.objects.FiberBundle))

        fbs = fastpli.objects.FiberBundles([[a]])
        self.assertTrue(isinstance(fbs, fastpli.objects.FiberBundles))
        fbs = fastpli.objects.FiberBundles(f)
        self.assertTrue(isinstance(fbs, fastpli.objects.FiberBundles))
        fbs = fastpli.objects.FiberBundles([f, f])
        self.assertTrue(isinstance(fbs, fastpli.objects.FiberBundles))
        fbs = fastpli.objects.FiberBundles(fbs)
        self.assertTrue(isinstance(fbs, fastpli.objects.FiberBundles))

        fb = fastpli.objects.FiberBundle([[[0, 0, 0, 1], [1, 1, 1, 1],
                                           [2, 2, 2, 1]],
                                          [[1, 0, 0, 1], [1, 1, 1, 1],
                                           [2, 2, 2, 1]]])
        for f in fb:
            self.assertTrue(isinstance(f, fastpli.objects.Fiber))
            self.assertTrue(isinstance(f._data, np.ndarray))

        fbs = fastpli.objects.FiberBundles([[[[0, 0, 0, 1], [1, 1, 1, 1],
                                              [2, 2, 2, 1]],
                                             [[1, 0, 0, 1], [1, 1, 1, 1],
                                              [2, 2, 2, 1]]],
                                            [[[0, 1, 2, 3], [1, 2, 3, 4],
                                              [2, 4, 5, 5]],
                                             [[1, 1, 2, 3], [1, 2, 3, 4],
                                              [2, 4, 5, 5]],
                                             [[1, 1, 2, 3], [1, 2, 3, 4],
                                              [2, 4, 5, 5]]]])

        for fb in fbs:
            self.assertTrue(isinstance(fb, fastpli.objects.FiberBundle))
            for f in fb:
                self.assertTrue(isinstance(f, fastpli.objects.Fiber))
                self.assertTrue(isinstance(f._data, np.ndarray))

    def test_type(self):
        self.assertTrue(isinstance(self.fiber[:], np.ndarray))
        self.assertTrue(self.fiber[:].dtype == float)
        self.assertTrue(
            fastpli.objects.Fiber([[1, 1, 1, 1]], np.float32).dtype ==
            np.float32)

    def test_layers(self):
        fastpli.objects.FiberBundle(self.fiber_bundle,
                                    [(0.333, -0.004, 10, 'p'),
                                     (0.666, 0, 5, 'b'), (1.0, 0.004, 1, 'r')])
        fastpli.objects.FiberBundles(self.fiber_bundles,
                                     [[(0.333, -0.004, 10, 'p'),
                                       (0.666, 0, 5, 'b'),
                                       (1.0, 0.004, 1, 'r')]])

        fb = fastpli.objects.FiberBundle([[[0, 0, 0, 1], [1, 1, 1, 1],
                                           [2, 2, 2, 1]],
                                          [[1, 0, 0, 1], [1, 1, 1, 1],
                                           [2, 2, 2, 1]]])
        fb = fastpli.objects.FiberBundle(fb, [(0.333, -0.004, 10, 'p'),
                                              (0.666, 0, 5, 'b'),
                                              (1.0, 0.004, 1, 'r')])

        fbs = [[[[0, 0, 0, 1], [1, 1, 1, 1], [2, 2, 2, 1]],
                [[1, 0, 0, 1], [1, 1, 1, 1], [2, 2, 2, 1]]],
               [[[0, 1, 2, 3], [1, 2, 3, 4], [2, 4, 5, 5]],
                [[1, 1, 2, 3], [1, 2, 3, 4], [2, 4, 5, 5]],
                [[1, 1, 2, 3], [1, 2, 3, 4], [2, 4, 5, 5]]]]
        fbs = fastpli.objects.FiberBundles(fbs,
                                           [[(0.333, -0.004, 10, 'p'),
                                             (0.666, 0, 5, 'b'),
                                             (1.0, 0.004, 1, 'r')]] * len(fbs))

    def test_resize(self):
        fiber = self.fiber.scale(10)
        self.assertTrue(np.array_equal(fiber[:], self.fiber[:] * 10))

        fb = self.fiber_bundle.scale(10)
        for f in fb:
            self.assertTrue(np.array_equal(f[:], self.fiber[:] * 10))

        fbs = self.fiber_bundles.scale(10)
        for fb in fbs:
            for f in fb:
                self.assertTrue(np.array_equal(f[:], self.fiber[:] * 10))

        fiber = self.fiber.scale(10, mode='points')
        self.assertTrue(np.array_equal(fiber[:, :-2], self.fiber[:, :-2] * 10))
        self.assertTrue(np.array_equal(fiber[:, -1], self.fiber[:, -1]))

        fiber = self.fiber.scale(10, mode='radii')
        self.assertTrue(np.array_equal(fiber[:, :-2], self.fiber[:, :-2]))
        self.assertTrue(np.array_equal(fiber[:, -1], self.fiber[:, -1] * 10))

    def test_rotation(self):

        fiber = self.fiber.rotate(fastpli.tools.rotation.x(0))
        self.assertTrue(np.array_equal(self.fiber[:], fiber[:]))

        fiber = self.fiber.rotate(fastpli.tools.rotation.x(np.deg2rad(90)))

        self.assertTrue(
            np.allclose(fiber[:], np.array([[0, 0, 0, 1], [1, -1, 1, 2]])))

        fiber = self.fiber.rotate(fastpli.tools.rotation.x(np.deg2rad(90)),
                                  [1, 1, 1])
        self.assertTrue(
            np.allclose(fiber[:], np.array([[0, 2, 0, 1], [1, 1, 1, 2]])))

        fiber_bundle = self.fiber_bundle.rotate(
            fastpli.tools.rotation.x(np.deg2rad(90)), [1, 1, 1])

        self.assertTrue(len(fiber_bundle) == len(self.fiber_bundle))
        for f in fiber_bundle:
            self.assertTrue(
                np.allclose(f[:], np.array([[0, 2, 0, 1], [1, 1, 1, 2]])))

        for fb in self.fiber_bundles:
            for f in fb:
                fiber = f.rotate(fastpli.tools.rotation.x(np.deg2rad(90)),
                                 [1, 1, 1])
                self.assertTrue(
                    np.allclose(fiber[:], np.array([[0, 2, 0, 1], [1, 1, 1,
                                                                   2]])))

    def test_translate(self):
        fiber = self.fiber.translate([1, 1, 1])
        self.assertTrue(
            np.array_equal(fiber[:, :3],
                           self.fiber[:, :3] + np.array([1, 1, 1])))
        self.assertTrue(np.array_equal(fiber[:, -1], self.fiber[:, -1]))

        fiber_bundle = self.fiber_bundle.translate([1, 1, 1])

        self.assertTrue(len(fiber_bundle) == len(self.fiber_bundle))
        for f in fiber_bundle:
            self.assertTrue(
                np.array_equal(fiber[:, :3],
                               self.fiber[:, :3] + np.array([1, 1, 1])))
            self.assertTrue(np.array_equal(f[:, -1], self.fiber[:, -1]))

        for fb in self.fiber_bundles:
            for f in fb:
                fiber = f.translate([1, 1, 1])
                self.assertTrue(
                    np.array_equal(fiber[:, :3],
                                   self.fiber[:, :3] + np.array([1, 1, 1])))
                self.assertTrue(np.array_equal(f[:, -1], self.fiber[:, -1]))

    def test_apply(self):
        # Fiber
        fiber = fastpli.objects.Fiber([[0, 0, 0, 1], [1, 1, 1, 2]], dtype=float)
        fiber_ = fiber.apply(lambda x: x + 1)
        self.assertTrue(isinstance(fiber_, fastpli.objects.Fiber))
        self.assertTrue(np.array_equal(fiber[:] + 1, fiber_[:]))

        fiber_ = fiber.apply_to_points(lambda x: x + 1)
        self.assertTrue(isinstance(fiber_, fastpli.objects.Fiber))
        self.assertTrue(np.array_equal(fiber[:, :-1] + 1, fiber_[:, :-1]))
        self.assertTrue(np.array_equal(fiber[:, -1], fiber_[:, -1]))

        fiber_ = fiber.apply_to_radii(lambda x: x + 1)
        self.assertTrue(isinstance(fiber_, fastpli.objects.Fiber))
        self.assertTrue(np.array_equal(fiber[:, :-1], fiber_[:, :-1]))
        self.assertTrue(np.array_equal(fiber[:, -1] + 1, fiber_[:, -1]))

        # FiberBundle
        fb = fastpli.objects.FiberBundle([[0, 0, 0, 1], [1, 1, 1, 2]],
                                         dtype=float)
        fb_ = fb.apply(lambda x: x + 1)
        self.assertTrue(isinstance(fb_, fastpli.objects.FiberBundle))
        self.assertTrue(np.array_equal(fb[0][:] + 1, fb_[0][:]))

        fb_ = fb.apply_to_points(lambda x: x + 1)
        self.assertTrue(isinstance(fb_, fastpli.objects.FiberBundle))
        self.assertTrue(np.array_equal(fb[0][:, :-1] + 1, fb_[0][:, :-1]))
        self.assertTrue(np.array_equal(fb[0][:, -1], fb_[0][:, -1]))

        fb_ = fb.apply_to_radii(lambda x: x + 1)
        self.assertTrue(isinstance(fb_, fastpli.objects.FiberBundle))
        self.assertTrue(np.array_equal(fb[0][:, :-1], fb_[0][:, :-1]))
        self.assertTrue(np.array_equal(fb[0][:, -1] + 1, fb_[0][:, -1]))

        # FiberBundles
        fbs = fastpli.objects.FiberBundles([[[0, 0, 0, 1], [1, 1, 1, 2]]],
                                           dtype=float)
        fbs_ = fbs.apply(lambda x: x + 1)

        self.assertTrue(isinstance(fbs_, fastpli.objects.FiberBundles))
        self.assertTrue(np.array_equal(fbs[0][0][:] + 1, fbs_[0][0][:]))

        fbs_ = fbs.apply_to_points(lambda x: x + 1)
        self.assertTrue(isinstance(fbs_, fastpli.objects.FiberBundles))
        self.assertTrue(
            np.array_equal(fbs[0][0][::, :-1] + 1, fbs_[0][0][:, :-1]))
        self.assertTrue(np.array_equal(fbs[0][0][:, -1], fbs_[0][0][:, -1]))

        fbs_ = fbs.apply_to_radii(lambda x: x + 1)
        self.assertTrue(isinstance(fbs_, fastpli.objects.FiberBundles))
        self.assertTrue(np.array_equal(fbs[0][0][:, :-1], fbs_[0][0][:, :-1]))
        self.assertTrue(np.array_equal(fbs[0][0][:, -1] + 1, fbs_[0][0][:, -1]))

    def test_cut(self):
        fiber = fastpli.objects.Fiber([[0, 0, 0, 1], [1, 1, 1, 2]], dtype=float)
        fibers = fiber.cut([[-10] * 3, [10] * 3])
        self.assertTrue(len(fibers) == 1)
        self.assertTrue(np.array_equal(fibers[0][:], fiber[:]))

        fiber = fastpli.objects.Fiber([[0, 0, 0, 1], [10, 10, 10, 2]])
        fibers = fiber.cut([[-5] * 3, [5] * 3])
        self.assertTrue(len(fibers) == 1)
        self.assertTrue(np.array_equal(fibers[0][:], fiber[:]))

        fiber = fastpli.objects.Fiber([[0, 0, 0, 1], [10, 10, 10, 2],
                                       [100, 100, 100, 2]])
        fibers = fiber.cut([[-5] * 3, [5] * 3])
        self.assertTrue(len(fibers) == 1)
        self.assertTrue(fibers[0].shape[0] == 2)
        self.assertTrue(not np.array_equal(fibers[0][:], fiber[:]))

        fiber = fastpli.objects.Fiber([[0, 0, 0, 1], [10, 10, 10, 2],
                                       [100, 100, 100, 2], [10, 10, 10, 2],
                                       [0, 0, 0, 1]])
        fibers = fiber.cut([[-5] * 3, [5] * 3])
        self.assertTrue(len(fibers) == 2)
        self.assertTrue(fibers[0].shape[0] == 2)
        self.assertTrue(fibers[1].shape[0] == 2)
        self.assertTrue(not np.array_equal(fibers[0][:], fiber[:]))
        self.assertTrue(not np.array_equal(fibers[1][:], fiber[:]))

        fiber_bundle = fastpli.objects.FiberBundle(fiber)
        cut_fb = fiber_bundle.cut([[-5] * 3, [5] * 3])
        fibers = cut_fb
        self.assertTrue(len(fibers) == 2)
        self.assertTrue(fibers[0].shape[0] == 2)
        self.assertTrue(fibers[1].shape[0] == 2)
        self.assertTrue(not np.array_equal(fibers[0][:], fiber[:]))
        self.assertTrue(not np.array_equal(fibers[1][:], fiber[:]))

        fiber_bundles = fastpli.objects.FiberBundles(fiber)
        cut_fbs = fiber_bundles.cut([[-5] * 3, [5] * 3])
        fibers = cut_fbs[0]
        self.assertTrue(len(cut_fbs) == 1)
        self.assertTrue(len(fibers) == 2)
        self.assertTrue(fibers[0].shape[0] == 2)
        self.assertTrue(fibers[1].shape[0] == 2)
        self.assertTrue(not np.array_equal(fibers[0][:], fiber[:]))
        self.assertTrue(not np.array_equal(fibers[1][:], fiber[:]))

        fiber = fastpli.objects.Fiber([[0, 0, 0, 1], [10, 10, 10, 2]])
        fibers = fiber.cut([[5] * 3, [6] * 3])
        self.assertTrue(np.array_equal(fibers[0][:], fiber[:]))


if __name__ == '__main__':
    unittest.main()
