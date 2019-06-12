from distutils.core import setup

import sys
import setuptools

if sys.version_info < (3, 0):
    sys.exit('Sorry, Python < 3.0 is not supported')

setup(
    name='fastpli',
    version='@GIT_DESCRIBE_LOG@-${CMAKE_BUILD_TYPE}',
    description='Fiber Architecture Simulation Toolbox for PLI',
    long_description='',
    author='Felix Matuschke',
    author_email='f.matuschke@fz-juelich.de',
    url=
    'http://www.fz-juelich.de/inm/inm-1/EN/Forschung/Fibre%20Architecture/Fibre%20Architecture_node.html',
    install_requires=[
        'numpy', 'numba', 'pymp-pypi', 'mpi4py', 'scipy', 'pillow', 'h5py'
    ],
    # test_suite='${CMAKE_SOURCE_DIR}/tests', # use make test instead
    zip_safe=False,
    packages=setuptools.find_packages(),
    package_data={'': ['*.so']})
