import unittest
import numpy as np

import fastpli.model.sandbox as sb


class MainTest(unittest.TestCase):

    def setUp(self):
        pass

    def test_triangular_grid(self):
        seeds = sb.seeds.triangular_grid(0, 0, 1, endpoint=False)
        self.assertTrue(seeds.size == 0)
        seeds = sb.seeds.triangular_grid(2, 0, 1, endpoint=False)
        self.assertTrue(seeds.size == 0)

        seeds = sb.seeds.triangular_grid(2, 0, 1, endpoint=True)
        self.assertTrue(seeds.shape[0] == 3)
        self.assertTrue(np.all(seeds[:, 1] == 0))
        self.assertTrue(np.all(seeds[:, 0] == np.array([0, 1, 2])))

        seeds = sb.seeds.triangular_grid(2, 2, 1, endpoint=True)
        self.assertTrue(seeds.shape[0] == 8)

    def test_crop_rectangle(self):
        seeds = sb.seeds.triangular_grid(2, 2, 1)
        new_seeds = sb.seeds.crop_rectangle(2, 2, seeds)
        self.assertTrue(np.array_equal(seeds, new_seeds))

        seeds = sb.seeds.triangular_grid(0, 0, 1)
        new_seeds = sb.seeds.crop_rectangle(2, 2, seeds, 1)
        self.assertTrue(new_seeds.size == 0)

        new_seeds = sb.seeds.crop_rectangle([-1, -1], [1, 1], seeds, 0)
        self.assertTrue(np.array_equal(new_seeds, [[0, 0]]))

        seeds = sb.seeds.triangular_grid(2, 2, 1)
        new_seeds = sb.seeds.crop_rectangle(2,
                                            2,
                                            seeds,
                                            radii=[1] * seeds.shape[0])
        self.assertTrue(new_seeds.size < seeds.size)

    def test_crop_circle(self):
        seeds = sb.seeds.triangular_grid(2, 2, 1)
        new_seeds = sb.seeds.crop_circle(100, seeds)
        self.assertTrue(np.array_equal(seeds, new_seeds))

        new_seeds = sb.seeds.crop_circle(1, seeds)
        self.assertTrue(
            np.array_equal(new_seeds,
                           [[0, 0], [1, 0], [0.5, np.sqrt(3) / 2]]))

        new_seeds = sb.seeds.crop_circle(np.sqrt(2), seeds, [1, 1])
        self.assertTrue(np.array_equal(seeds, new_seeds))

        new_seeds = sb.seeds.crop_circle(1, seeds, center=[0, 0], radii=1)
        self.assertTrue(np.array_equal(new_seeds, [[0, 0]]))

        new_seeds = sb.seeds.crop_circle(1,
                                         seeds,
                                         center=[0, 0],
                                         radii=[1] * seeds.shape[0])
        self.assertTrue(np.array_equal(new_seeds, [[0, 0]]))


if __name__ == '__main__':
    unittest.main()
