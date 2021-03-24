import setuptools
import sys

import src.fastpli.__version

version = '@GIT_DESCRIBE_LOG@'
i = version.find('-')
if i > 0:
    version = version[:i] + '.dev' + version[i + 1:] + '-${CMAKE_BUILD_TYPE}'

compiled_version = src.fastpli.__version.__libraries__.split(";")[0].split(
    "v")[-1]
sys_version = f"{sys.version_info[0]}.{sys.version_info[1]}"

if not compiled_version.startswith(sys_version):
    raise SystemError(
        f"Unsuported compiled python library version: compiled: {compiled_version} != system: {sys_version}. Specify Python version in cmake configure file"
    )

setuptools.setup(
    name='fastpli',
    version=version,
    description='Fiber Architecture Simulation Toolbox for PLI',
    long_description='',
    author='Felix Matuschke',
    author_email='f.matuschke@fz-juelich.de',
    url='http://www.fz-juelich.de/inm/inm-1/EN/Forschung/\
        Fibre%20Architecture/Fibre%20Architecture_node.html',
    python_requires='>3.6.0',
    install_requires=[
        'numpy>=1.19', 'numba>=0.52', 'scipy>=1.5.4', 'h5py>=3.1'
    ],
    zip_safe=False,
    packages=setuptools.find_packages('src'),
    package_dir={'': 'src'},
    package_data={'': ['*.so']},
)
