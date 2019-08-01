import unittest


class MainTest(unittest.TestCase):

    def test_import(self):
        import fastpli
        import fastpli.analysis
        import fastpli.io
        import fastpli.model
        import fastpli.objects
        import fastpli.simulation
        import fastpli.tools

        print(fastpli.__version__)


if __name__ == '__main__':
    unittest.main()
