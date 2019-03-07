import os
import unittest
import numpy as np

import fastpli


class MainTest(unittest.TestCase):

    def setUp(self):
        self._test_fiber = fastpli.objects.Fiber([0, 0, 0, 0, 0, 1], [1, 2])
        self._test_fiberbundles = [[self._test_fiber]]
        self.solver = fastpli.model.Solver()
        self.solver.fiber_bundles = self._test_fiberbundles

    def test_number_of_fibers(self):
        self.solver.set_parameters(drag=0, obj_min_radius=0, obj_mean_length=0)
        self.solver.step()
        fb = self.solver.fiber_bundles
        self.assertTrue(
            np.array_equal(self._test_fiber.points, fb[0][0].points))
        self.assertTrue(np.array_equal(self._test_fiber.radii, fb[0][0].radii))

    def test_fiber_bundle_property(self):
        self.solver.parameters = (0, 0, 0)
        _ = self.solver.parameters
        self.solver.fiber_bundles = [
            [fastpli.objects.Fiber([0, 0, 0, 0, 0, 2], [1, 3])]]
        test = self.solver.fiber_bundles

    def test_split(self):
        self.solver.set_parameters(
            drag=0, obj_min_radius=0, obj_mean_length=0.5)
        self.solver.step()
        fb = self.solver.fiber_bundles
        self.assertTrue(fb[0][0].radii[1] == 1.5)

    def test_combine(self):
        self.fiber = fastpli.objects.Fiber(
            [0, 0, 0, 0, 0, 1, 0, 0, 2], [1, 1, 1])
        self._test_fiberbundles = [[self._test_fiber]]
        self.solver.fiber_bundles = self._test_fiberbundles
        self.solver.set_parameters(drag=0, obj_min_radius=0, obj_mean_length=2)
        self.solver.step()
        fb = self.solver.fiber_bundles
        self.assertTrue(fb[0][0].radii.shape[0] == 2)

        self.solver.set_parameters(
            drag=0,
            obj_min_radius=0,
            obj_mean_length=20)
        self.solver.step()
        fb = self.solver.fiber_bundles
        self.assertTrue(fb[0][0].radii.shape[0] == 2)

    def test_fibers(self):
        self.solver.set_parameters(drag=0, obj_min_radius=0, obj_mean_length=2)

        fiber_0 = fastpli.objects.Fiber([0, 0, 0, 0, 0, 1], [1, 2])
        fiber_1 = fastpli.objects.Fiber([0, 0, 0.1, 0, 0, 1.1], [1, 2])
        self.solver.fiber_bundles = [[fiber_0, fiber_1]]
        self.solver.step()
        fbs = self.solver.fiber_bundles

        self.assertFalse(np.array_equal(fiber_0.points, fbs[0][0].points))
        self.assertTrue(np.array_equal(fiber_0.radii, fbs[0][0].radii))
        self.assertFalse(np.array_equal(fiber_1.points, fbs[0][1].points))
        self.assertTrue(np.array_equal(fiber_1.radii, fbs[0][1].radii))

    def test_fiber_bundles(self):
        self.solver.set_parameters(drag=0, obj_min_radius=0, obj_mean_length=2)

        fiber_0 = fastpli.objects.Fiber([0, 0, 0, 0, 0, 1], [1, 2])
        fiber_1 = fastpli.objects.Fiber([0, 0, 0.1, 0, 0, 1.1], [1, 2])
        self.solver.fiber_bundles = [[fiber_0], [fiber_1]]
        self.solver.step()
        fbs = self.solver.fiber_bundles

        self.assertFalse(np.array_equal(fiber_0.points, fbs[0][0].points))
        self.assertTrue(np.array_equal(fiber_0.radii, fbs[0][0].radii))
        self.assertFalse(np.array_equal(fiber_1.points, fbs[1][0].points))
        self.assertTrue(np.array_equal(fiber_1.radii, fbs[1][0].radii))

    def test_col_voi(self):
        self.solver.set_parameters(drag=0, obj_min_radius=0, obj_mean_length=2)

        fiber_0 = fastpli.objects.Fiber([0, 0, 0, 0, 0, 1], [1, 2])
        fiber_1 = fastpli.objects.Fiber([0, 0, 0.1, 0, 0, 1.1], [1, 2])
        self.solver.fiber_bundles = [[fiber_0], [fiber_1]]

        self.solver.col_voi = ([-10, -10, -10], [-9, -9, -9])
        self.solver.step()
        fbs = self.solver.fiber_bundles
        self.assertTrue(np.array_equal(fiber_0.points, fbs[0][0].points))
        self.assertTrue(np.array_equal(fiber_0.radii, fbs[0][0].radii))

        self.solver.col_voi = ([0, 0, 0], [1, 1, 1])
        self.solver.step()
        fbs = self.solver.fiber_bundles
        self.assertFalse(np.array_equal(fiber_0.points, fbs[0][0].points))
        self.assertTrue(np.array_equal(fiber_0.radii, fbs[0][0].radii))

    def test_openmp(self):
        self.solver.omp_num_threads = 2
        self.assertTrue(self.solver.omp_num_threads >= 0)

    def test_opengl(self):
        display = ""
        try:
            display = os.environ['DISPLAY']
        except BaseException:
            print("test_opengl: no display detected")
        if display:
            self.solver.draw_scene()
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
