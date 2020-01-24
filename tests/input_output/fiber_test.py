import unittest
import numpy as np
import h5py
import os

from fastpli.io import fiber
from fastpli.model.solver import Solver


class MainTest(unittest.TestCase):

    def setUp(self):
        self.fiber_bundles = [[[[0, 0, 0, 1], [1, 1, 1, 1], [2, 2, 2, 1]],
                               [[1, 0, 0, 1], [1, 1, 1, 1], [2, 2, 2, 1]]],
                              [[[0, 1, 2, 3], [1, 2, 3, 4], [2, 4, 5, 5]],
                               [[1, 1, 2, 3], [1, 2, 3, 4], [2, 4, 5, 5]]]]

        self.solver = Solver()
        self.solver.fiber_bundles = self.fiber_bundles

    def test_h5(self):
        with h5py.File('/tmp/fastpli.test.h5', 'w') as h5f:
            fiber.save_h5(h5f, self.solver.fiber_bundles)

        with h5py.File('/tmp/fastpli.test.h5', 'r') as h5f:
            fbs = fiber.load_h5(h5f)
            for fb_a, fb_b in zip(fbs, self.solver.fiber_bundles):
                for f_a, f_b in zip(fb_a, fb_b):
                    self.assertTrue(np.alltrue(f_a == f_b))

        with h5py.File('/tmp/fastpli.test.h5', 'w') as h5f:
            h5g = h5f.create_group('test')
            fiber.save_h5(h5g, self.solver.fiber_bundles)

        with h5py.File('/tmp/fastpli.test.h5', 'r') as h5f:
            fbs = fiber.load_h5(h5f['test'])
            for fb_a, fb_b in zip(fbs, self.solver.fiber_bundles):
                for f_a, f_b in zip(fb_a, fb_b):
                    self.assertTrue(np.alltrue(f_a == f_b))

        self.addCleanup(os.remove, '/tmp/fastpli.test.h5')

    def test_dat(self):
        with open('/tmp/fastpli.test.dat', 'w') as file:
            fiber.save_dat(file, self.solver.fiber_bundles)

        with open('/tmp/fastpli.test.dat', 'r') as file:
            fbs = fiber.load_dat(file)
            for fb_a, fb_b in zip(fbs, self.solver.fiber_bundles):
                for f_a, f_b in zip(fb_a, fb_b):
                    self.assertTrue(np.alltrue(f_a == f_b))

        self.addCleanup(os.remove, '/tmp/fastpli.test.dat')

    def test_h5_dat(self):
        with h5py.File('/tmp/fastpli.test.h5', 'w') as h5f:
            fiber.save_h5(h5f, self.solver.fiber_bundles)

        with open('/tmp/fastpli.test.dat', 'w') as file:
            fiber.save_dat(file, self.solver.fiber_bundles)

        with h5py.File('/tmp/fastpli.test.h5', 'r') as h5f:
            with open('/tmp/fastpli.test.dat', 'r') as file:
                fbs_h5 = fiber.load_h5(h5f)
                fbs_dat = fiber.load_dat(file)
                for fb_a, fb_b in zip(fbs_h5, fbs_dat):
                    for f_a, f_b in zip(fb_a, fb_b):
                        self.assertTrue(np.alltrue(f_a == f_b))

        self.addCleanup(os.remove, '/tmp/fastpli.test.h5')
        self.addCleanup(os.remove, '/tmp/fastpli.test.dat')

    def test_load_save(self):
        fiber.save('/tmp/fastpli.test.dat', self.solver.fiber_bundles,
                   'fiber_bundles'
                   'w-')
        fiber.save('/tmp/fastpli.test.h5', self.solver.fiber_bundles,
                   'fiber_bundles', 'a')

        fbs_h5 = fiber.load('/tmp/fastpli.test.h5', 'fiber_bundles')
        fbs_dat = fiber.load('/tmp/fastpli.test.dat', 'fiber_bundles')
        for fb_a, fb_b in zip(fbs_h5, fbs_dat):
            for f_a, f_b in zip(fb_a, fb_b):
                self.assertTrue(np.alltrue(f_a == f_b))

        self.addCleanup(os.remove, '/tmp/fastpli.test.h5')
        self.addCleanup(os.remove, '/tmp/fastpli.test.dat')


if __name__ == '__main__':
    unittest.main()
