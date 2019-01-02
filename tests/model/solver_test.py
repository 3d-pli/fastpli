import unittest
import numpy as np

import fastpli


class MainTest(unittest.TestCase):

    def setUp(self):
        self.fiber = fastpli.objects.Fiber([0, 0, 0, 0, 0, 1], [1, 2])
        self.fiberbundles = [[self.fiber]]
        self.solver = fastpli.model.Solver()
        self.solver.set_fiber_bundles(self.fiberbundles)

    def test_number_of_fibers(self):
        self.solver.set_parameters(drag=0, obj_min_radius=0, obj_mean_length=0)
        self.solver.step()
        fb = self.solver.get_fiber_bundles()
        self.assertTrue(np.array_equal(self.fiber.points, fb[0][0].points))
        self.assertTrue(np.array_equal(self.fiber.radii, fb[0][0].radii))

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
        fb = self.solver.get_fiber_bundles()
        self.assertTrue(fb[0][0].radii[1] == 1.5)

    def test_combine(self):
        self.fiber = fastpli.objects.Fiber(
            [0, 0, 0, 0, 0, 1, 0, 0, 2], [1, 1, 1])
        self.fiberbundles = [[self.fiber]]
        self.solver.set_fiber_bundles(self.fiberbundles)
        self.solver.set_parameters(drag=0, obj_min_radius=0, obj_mean_length=2)
        self.solver.step()
        fb = self.solver.get_fiber_bundles()
        self.assertTrue(fb[0][0].radii.shape[0] == 2)

        self.solver.set_parameters(
            drag=0,
            obj_min_radius=0,
            obj_mean_length=20)
        self.solver.step()
        fb = self.solver.get_fiber_bundles()
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

        self.solver.set_col_voi([-10, -10, -10], [-9, -9, -9])
        self.solver.step()
        fbs = self.solver.fiber_bundles
        self.assertTrue(np.array_equal(fiber_0.points, fbs[0][0].points))
        self.assertTrue(np.array_equal(fiber_0.radii, fbs[0][0].radii))

        self.solver.set_col_voi([0, 0, 0], [1, 1, 1])
        self.solver.step()
        fbs = self.solver.fiber_bundles
        self.assertFalse(np.array_equal(fiber_0.points, fbs[0][0].points))
        self.assertTrue(np.array_equal(fiber_0.radii, fbs[0][0].radii))

    def test_openmp(self):
        i = self.solver.set_omp_num_threads(2)
        self.assertTrue(i >= 0)


if __name__ == '__main__':
    unittest.main()
