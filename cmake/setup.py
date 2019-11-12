import os
import sys
import setuptools

if sys.version_info < (3, 6):
    sys.exit('Python >= 3.6 required')

# copy src files into build folder
os.system(
    'cp -al --remove-destination ${CMAKE_SOURCE_DIR}/src/fastpli ${CMAKE_CURRENT_BINARY_DIR}/'
)

setuptools.setup(
    name='fastpli',
    version='@GIT_DESCRIBE_LOG@-${CMAKE_BUILD_TYPE}',
    description='Fiber Architecture Simulation Toolbox for PLI',
    long_description='',
    author='Felix Matuschke',
    author_email='f.matuschke@fz-juelich.de',
    url=
    'http://www.fz-juelich.de/inm/inm-1/EN/Forschung/Fibre%20Architecture/Fibre%20Architecture_node.html',
    install_requires=['numpy', 'numba', 'mpi4py', 'scipy', 'h5py'],
    # test_suite='' # use make test instead
    zip_safe=False,
    package_dir={'': os.path.relpath('${CMAKE_CURRENT_BINARY_DIR}')},
    packages=setuptools.find_packages(
        os.path.relpath('${CMAKE_CURRENT_BINARY_DIR}')),
    package_data={'': ['*.so']})
