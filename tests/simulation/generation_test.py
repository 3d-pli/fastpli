import unittest
import numpy as np

import fastpli
from fastpli.simulation.helper import TupleList2Layer


class MainTest(unittest.TestCase):

    def setUp(self):
        self.points = [0, 0, 0, 1, 1, 1, 2, 2, 2]
        self.radii = [1, 1, 1]
        self.fiber_bundles = [[fastpli.objects.Fiber(self.points, self.radii)]]
        fiber_prop = [[(0.333, 0.004, 10, 'p'), (
    0.666, -0.004, 5, 'b'), (1.0, 0.004, 1, 'r')]]
        self.fiber_prop = TupleList2Layer(fiber_prop)

        self.generator = fastpli.simulation.generation.Generator()
        self.generator.set_volume([10, 10, 10], [10, 10, 10], [0, 0, 0], 0.2)
        self.generator.set_fiber_bundles(self.fiber_bundles, self.fiber_prop)

    def test_return_dimension(self):

        fiber_bundles = [[fastpli.objects.Fiber([1, 3, 0, 1, 3, 7], [2, 2])]]
        fiber_prop = [[(1, 0, 0, 'p')]]
        fiber_prop = TupleList2Layer(fiber_prop)
        self.generator.set_fiber_bundles(fiber_bundles, fiber_prop)
        self.generator.set_volume([3, 5, 7], [3, 5, 7], [0, 0, 0], 1)

        label_field, vec_field, tissue_properties = self.generator.run_generation(
            )

        self.assertTrue(np.array_equal(label_field.shape, [3, 5, 7]))
        self.assertTrue(np.array_equal(vec_field.shape, [3, 5, 7, 3]))
        self.assertTrue(
            np.array_equal(label_field.shape, vec_field.shape[0:3]))
        self.assertTrue(np.array_equal(label_field, np.array(vec_field[:, :,:, 2], dtype=np.int16)))

    def test_generator(self):
        label_field_0, vec_field, tissue_properties = self.generator.run_generation(
            only_label=True)
        self.assertTrue(vec_field.size == 0)

        label_field_1, vec_field, tissue_properties = self.generator.run_generation(
            only_label=False)
        self.assertTrue(np.array_equal(label_field_0, label_field_1))
        self.assertTrue(
            np.array_equal(label_field_1.shape, vec_field.shape[0:3]))


if __name__ == '__main__':
    unittest.main()
