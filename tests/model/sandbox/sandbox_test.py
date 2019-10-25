import unittest
import numpy as np

import fastpli.model.sandbox as sb


class MainTest(unittest.TestCase):

    def setUp(self):
        pass

    def test_seeds(self):
        seeds = sb.seeds.triangle_grid(0, 0, 1)
        self.assertTrue(np.all(seeds == np.array([0, 0])))

        seeds = sb.seeds.triangle_grid(2, 0, 1)
        self.assertTrue(seeds.shape[0] == 3)
        self.assertTrue(np.all(seeds[:, 1] == 0))
        self.assertTrue(np.all(seeds[:, 0] == np.array([0, 1, 2])))

        seeds = sb.seeds.triangle_grid(2, 2, 1)
        self.assertTrue(seeds.shape[0] == 8)

    def test_crop_rectangle(self):
        seeds = sb.seeds.triangle_grid(2, 2, 1)
        new_seeds = sb.seeds.crop_rectangle(2, 2, seeds)
        self.assertTrue(np.all(seeds == new_seeds))

        seeds = sb.seeds.triangle_grid(0, 0, 1)
        new_seeds = sb.seeds.crop_rectangle(2, 2, seeds, 1)
        self.assertTrue(new_seeds.size == 0)


if __name__ == '__main__':
    unittest.main()
