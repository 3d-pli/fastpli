import unittest
import numpy as np

import fastpli.model.solver


class MainTest(unittest.TestCase):

    def setUp(self):
        self._test_fiber = np.array([[0, 0, 0, 1], [0, 0, 1, 2]])
        self._test_fiberbundles = [[self._test_fiber]]
        self.solver = fastpli.model.solver.Solver()
        self.solver.fiber_bundles = self._test_fiberbundles

    def test_dict(self):
        self.solver.as_dict()

    def test_number_of_fibers(self):
        self.solver.drag = 0
        self.solver.obj_min_radius = 0
        self.solver.obj_mean_length = 0
        self.solver.step()
        fbs = self.solver.fiber_bundles
        self.assertTrue(fbs[0][0].shape[0] == 2)
        self.assertTrue(np.array_equal(self._test_fiber, fbs[0][0]))

    def test_set_fiber_bundle(self):
        self.solver.fiber_bundles = [[np.array([[0, 0, 0, 1], [0, 0, 2, 3]])]]
        _ = self.solver.fiber_bundles

    def test_split(self):
        self.solver.drag = 0
        self.solver.obj_min_radius = 0
        self.solver.obj_mean_length = 0.5
        self.solver.step()
        fbs = self.solver.fiber_bundles
        self.assertTrue(fbs[0][0].shape[0] == 2)
        self.assertTrue(fbs[0][0][1, -1] == 2)

        fbs = self.solver.apply_boundary_conditions(n_max=1)
        self.assertTrue(fbs[0][0].shape[0] == 3)
        self.assertTrue(np.isclose(fbs[0][0][1, -1], 1.5))

    def test_combine(self):
        self.fiber = np.array([[0, 0, 0, 1], [0, 0, 1, 1], [0, 0, 2, 1]])
        self._test_fiberbundles = [[self._test_fiber]]
        self.solver.fiber_bundles = self._test_fiberbundles
        self.solver.drag = 0
        self.solver.obj_min_radius = 0
        self.solver.obj_mean_length = 2
        self.solver.step()
        fbs = self.solver.fiber_bundles
        self.assertTrue(fbs[0][0].shape[0] == 2)

        self.solver.drag = 0
        self.solver.obj_min_radius = 0
        self.solver.obj_mean_length = 20
        self.solver.step()
        fbs = self.solver.fiber_bundles
        self.assertTrue(fbs[0][0].shape[0] == 2)

    def test_fibers(self):
        self.solver.drag = 0
        self.solver.obj_min_radius = 0
        self.solver.obj_mean_length = 2

        fiber_0 = np.array([[0, 0, 0, 1], [0, 0, 1, 2]])
        fiber_1 = np.array([[0, 0, 0.1, 1], [0, 0, 1.1, 2]])
        self.solver.fiber_bundles = [[fiber_0, fiber_1]]
        self.solver.step()
        fbs = self.solver.fiber_bundles

        self.assertFalse(np.array_equal(fiber_0[:, :3], fbs[0][0][:, :3]))
        self.assertTrue(np.array_equal(fiber_0[:, -1], fbs[0][0][:, -1]))
        self.assertFalse(np.array_equal(fiber_1[:, :3], fbs[0][1][:, :3]))
        self.assertTrue(np.array_equal(fiber_1[:, -1], fbs[0][1][:, -1]))

    def test_fiber_bundles(self):
        self.solver.drag = 0
        self.solver.obj_min_radius = 0
        self.solver.obj_mean_length = 2

        fiber_0 = np.array([[0, 0, 0, 1], [0, 0, 1, 2]])
        fiber_1 = np.array([[0, 0, 0.1, 1], [0, 0, 1.1, 2]])
        self.solver.fiber_bundles = [[fiber_0], [fiber_1]]
        self.solver.step()
        fbs = self.solver.fiber_bundles

        self.assertFalse(np.array_equal(fiber_0[:, :3], fbs[0][0][:, :3]))
        self.assertTrue(np.array_equal(fiber_0[:, -1], fbs[0][0][:, -1]))
        self.assertFalse(np.array_equal(fiber_1[:, :3], fbs[1][0][:, :3]))
        self.assertTrue(np.array_equal(fiber_1[:, -1], fbs[1][0][:, -1]))

    def test_col_voi(self):
        self.solver.drag = 0
        self.solver.obj_min_radius = 0
        self.solver.obj_mean_length = 2

        fiber_0 = np.array([[0, 0, 0, 1], [0, 0, 1, 2]])
        fiber_1 = np.array([[0, 0, 0.1, 1], [0, 0, 1.1, 2]])
        self.solver.fiber_bundles = [[fiber_0], [fiber_1]]

        self.solver.col_voi = ([-10, -10, -10], [-9, -9, -9])
        self.solver.step()
        fbs = self.solver.fiber_bundles
        self.assertTrue(np.array_equal(fiber_0, fbs[0][0]))
        self.assertTrue(np.array_equal(fiber_1, fbs[1][0]))

        self.solver.col_voi = ([0, 0, 0], [1, 1, 1])
        self.solver.step()
        fbs = self.solver.fiber_bundles
        self.assertFalse(np.array_equal(fiber_0[:, :3], fbs[0][0][:, :3]))
        self.assertTrue(np.array_equal(fiber_0[:, -1], fbs[0][0][:, -1]))

    def test_openmp(self):
        self.solver.omp_num_threads = 2
        self.assertTrue(self.solver.omp_num_threads >= 0)

    def test_opengl(self):
        self.solver.draw_scene()
        self.solver.draw_scene()
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
