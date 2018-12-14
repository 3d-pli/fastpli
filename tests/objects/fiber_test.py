import unittest
import fastpli
import numpy as np

class MainTest(unittest.TestCase):

    def setUp(self):
        self.points = np.array([0, 0, 0, 0, 0, 1, 0, 0, 2])
        self.radii = [1, 1, 1]
        self.fiber = fastpli.objects.Fiber(self.points, self.radii)

    def test_init(self):
        self.assertIsInstance(self.fiber.points, (np.ndarray))
        self.assertIsInstance(self.fiber.radii, (np.ndarray))
        self.assertTrue(np.array_equal(
            np.array(self.points, dtype=np.float32).reshape(-1, 3), self.fiber.points))
        self.assertTrue(
            np.array_equal(
                np.array(
                    self.radii,
                    dtype=np.float32),
                self.fiber.radii))

    # def test_rotation(self):
    #     fd.rotate_fiber()

if __name__ == '__main__':
    unittest.main()
