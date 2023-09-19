import setuptools
import sys

import src.fastpli.__version

# check python dev library version
compiled_version = src.fastpli.__version.__libraries__.split(";")[0].split(
    "v")[-1]
sys_version = f"{sys.version_info[0]}.{sys.version_info[1]}"

if not compiled_version.startswith(sys_version):
    raise SystemError(
        f"Unsuported compiled python library version: compiled: {compiled_version} != system: {sys_version}. Specify Python version in cmake configure file"
    )

setuptools.setup(
    name='fastpli',
    description='Fiber Architecture Simulation Toolbox for PLI',
    long_description='',
    author='Felix Matuschke',
    author_email='f.matuschke@fz-juelich.de',
    url='http://www.fz-juelich.de/inm/inm-1/EN/Forschung/\
        Fibre%20Architecture/Fibre%20Architecture_node.html',
    python_requires='>=3.10.0',
    install_requires=[
        'numpy>=1.24', 'numba>=0.57', 'scipy>=1.11.0', 'h5py>=3.9'
    ],
    zip_safe=False,
    packages=setuptools.find_packages('src'),
    package_dir={'': 'src'},
    package_data={'': ['*.so']},
    setuptools_git_versioning={
        "enabled": True,
    },
    setup_requires=["setuptools-git-versioning<2"],
)
