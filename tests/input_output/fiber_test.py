import unittest
import numpy as np
import h5py
import os
import tempfile

from fastpli.io import fiber_bundles
from fastpli.model.solver import Solver

TMP_FILE = os.path.join(tempfile.gettempdir(), "fastpli.test")


class MainTest(unittest.TestCase):
    def setUp(self):
        self.fiber_bundles = [[[[0, 0, 0, 1], [1, 1, 1, 1], [2, 2, 2, 1]],
                               [[1, 0, 0, 1], [1, 1, 1, 1], [2, 2, 2, 1]]],
                              [[[0, 1, 2, 3], [1, 2, 3, 4], [2, 4, 5, 5]],
                               [[1, 1, 2, 3], [1, 2, 3, 4], [2, 4, 5, 5]],
                               [[1, 1, 2, 3], [1, 2, 3, 4], [2, 4, 5, 5]]]]

        self.solver = Solver()
        self.solver.fiber_bundles = self.fiber_bundles

    def test_h5(self):
        with h5py.File(TMP_FILE + '.h5', 'w-') as h5f:
            fiber_bundles.save_h5(h5f, self.solver.fiber_bundles)

        with h5py.File(TMP_FILE + '.h5', 'r') as h5f:
            fbs = fiber_bundles.load_h5(h5f)
            for fb_a, fb_b in zip(fbs, self.solver.fiber_bundles):
                for f_a, f_b in zip(fb_a, fb_b):
                    self.assertTrue(np.alltrue(f_a == f_b))

        with h5py.File(TMP_FILE + '.h5', 'w') as h5f:
            h5g = h5f.create_group('test')
            fiber_bundles.save_h5(h5g, self.solver.fiber_bundles)

        with h5py.File(TMP_FILE + '.h5', 'r') as h5f:
            fbs = fiber_bundles.load_h5(h5f['test'])
            for fb_a, fb_b in zip(fbs, self.solver.fiber_bundles):
                for f_a, f_b in zip(fb_a, fb_b):
                    self.assertTrue(np.alltrue(f_a == f_b))

        self.addCleanup(os.remove, TMP_FILE + '.h5')

    def test_dat(self):
        with open(TMP_FILE + '.dat', 'w') as file:
            fiber_bundles.save_dat(file, self.solver.fiber_bundles)

        with open(TMP_FILE + '.dat', 'r') as file:
            fbs = fiber_bundles.load_dat(file)
            for fb_a, fb_b in zip(fbs, self.solver.fiber_bundles):
                for f_a, f_b in zip(fb_a, fb_b):
                    self.assertTrue(np.alltrue(f_a == f_b))

        self.addCleanup(os.remove, TMP_FILE + '.dat')

    def test_h5_dat(self):
        with h5py.File(TMP_FILE + '.h5', 'w-') as h5f:
            s = Solver()
            fiber_bundles.save_h5(h5f, s.fiber_bundles)

        with h5py.File(TMP_FILE + '.h5', 'w') as h5f:
            fiber_bundles.save_h5(h5f, self.solver.fiber_bundles)

        with open(TMP_FILE + '.dat', 'w') as file:
            fiber_bundles.save_dat(file, self.solver.fiber_bundles)

        with h5py.File(TMP_FILE + '.h5', 'r') as h5f:
            with open(TMP_FILE + '.dat', 'r') as file:
                fbs_h5 = fiber_bundles.load_h5(h5f)
                fbs_dat = fiber_bundles.load_dat(file)
                for fb_a, fb_b in zip(fbs_h5, fbs_dat):
                    for f_a, f_b in zip(fb_a, fb_b):
                        self.assertTrue(np.alltrue(f_a == f_b))

        self.addCleanup(os.remove, TMP_FILE + '.h5')
        self.addCleanup(os.remove, TMP_FILE + '.dat')

    def test_save_load(self):
        fiber_bundles.save(TMP_FILE + '.dat', self.solver.fiber_bundles,
                           'fiber_bundles'
                           'w-')
        fiber_bundles.save(TMP_FILE + '.h5', self.solver.fiber_bundles,
                           'fiber_bundles', 'a')

        fbs_h5 = fiber_bundles.load(TMP_FILE + '.h5', 'fiber_bundles')
        fbs_dat = fiber_bundles.load(TMP_FILE + '.dat', 'fiber_bundles')

        self.assertTrue(len(fbs_h5) == len(fbs_dat))

        for fb_a, fb_b in zip(fbs_h5, fbs_dat):
            self.assertTrue(len(fb_a) == len(fb_b))
            for f_a, f_b in zip(fb_a, fb_b):
                self.assertTrue(np.alltrue(f_a == f_b))

        self.addCleanup(os.remove, TMP_FILE + '.h5')
        self.addCleanup(os.remove, TMP_FILE + '.dat')


if __name__ == '__main__':
    unittest.main()
