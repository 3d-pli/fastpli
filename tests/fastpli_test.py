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
        print(fastpli.__git__name__)
        print(fastpli.__git__branch__)
        print(fastpli.__git__hash__)
        print(fastpli.__git__build__)
        print(fastpli.__compiler__)
        print(fastpli.__libraries__)


if __name__ == '__main__':
    unittest.main()
