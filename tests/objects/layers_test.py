import unittest

import fastpli.objects
import fastpli.tools


class MainTest(unittest.TestCase):

    def setUp(self):
        self.layer = (1, 0, 1, 'b')
        self.layers = [(0.5, 0, 1, 'b'), (0.75, 1, 1, 'r'), (1, -1, 1, 'p')]

    def test_init(self):
        fastpli.objects.Layer(self.layer)
        fastpli.objects.Layers(self.layers)


if __name__ == '__main__':
    unittest.main()
