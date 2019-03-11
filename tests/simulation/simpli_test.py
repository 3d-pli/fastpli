import unittest
import numpy as np

from fastpli.simulation import Simpli


class MainTest(unittest.TestCase):

    def setUp(self):
        self.fiber_bundles = [[[[0, 0, 0, 1], [1, 1, 1, 1], [2, 2, 2, 1]]]]
        self.fiber_bundles_properties = [[(0.333, 0.004, 10, 'p'),
                                          (0.666, -0.004, 5, 'b'),
                                          (1.0, 0.004, 1, 'r')]]

        self.simpli = Simpli()
        self.simpli.fiber_bundles = self.fiber_bundles
        self.simpli.fiber_bundles_properties = self.fiber_bundles_properties
        self.simpli.dim = [10, 10, 10]
        self.simpli.pixel_size = 0.2

    def test_return_dimension(self):
        self.simpli.fiber_bundles = [[[[1, 3, 0, 2], [1, 3, 7, 2]]]]
        self.simpli.fiber_bundles_properties = [[(1, 0, 0, 'p')]]
        self.simpli.dim = [3, 5, 7]
        self.simpli.origin = [0, 0, 0]
        self.simpli.pixel_size = 1

        label_field, vec_field, tissue_properties = self.simpli.GenerateTissue()

        self.assertTrue(np.array_equal(label_field.shape, [3, 5, 7]))
        self.assertTrue(np.array_equal(vec_field.shape, [3, 5, 7, 3]))
        self.assertTrue(np.array_equal(label_field.shape, vec_field.shape[0:3]))
        self.assertTrue(
            np.array_equal(label_field,
                           np.array(vec_field[:, :, :, 2], dtype=np.int16)))

    def test_generator(self):
        label_field_0, vec_field, tissue_properties = self.simpli.GenerateTissue(
            only_label=True)
        self.assertTrue(vec_field.size == 0)

        label_field_1, vec_field, tissue_properties = self.simpli.GenerateTissue(
            only_label=False)
        self.assertTrue(np.array_equal(label_field_0, label_field_1))
        self.assertTrue(
            np.array_equal(label_field_1.shape, vec_field.shape[0:3]))


if __name__ == '__main__':
    unittest.main()
