import unittest


class MainTest(unittest.TestCase):

    def test_import(self):
        import fastpli
        import fastpli.tools
        import fastpli.objects
        import fastpli.model
        import fastpli.simulation

        print(fastpli.__version__)


if __name__ == '__main__':
    unittest.main()
