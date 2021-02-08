<!--
________             ___________________________
___  __/_____ _________  /___  __ \__  /____  _/
__  /_ _  __ `/_  ___/  __/_  /_/ /_  /  __  /
_  __/ / /_/ /_(__  )/ /_ _  ____/_  /____/ /
/_/    \__,_/ /____/ \__/ /_/     /_____/___/
-->

# Fiber Architecture Simulation Toolbox for 3D-PLI

![fastpli-logo](logo.svg)

[fastPLI](https://github.com/3d-pli/fastpli) is an open source toolbox for modeling nerve fibers, simulating them in a [3D-PLI](https://dx.doi.org/10.3389%2Ffninf.2011.00034) microscope and the signal processing developed by the [fiber architecture group](https://www.fz-juelich.de/inm/inm-1/EN/Forschung/Fibre%20Architecture/Fibre%20Architecture_node.html) at the [Forschungszentrum Jülich](https://www.fz-juelich.de) - [INM1](https://www.fz-juelich.de/inm/inm-1/EN/Home/home_node.html).

It consists of three consecutive parts:

- [Sandbox](https://github.com/3d-pli/fastpli/wiki/Sandbox): building nerve fiber models.
- [Solver](https://github.com/3d-pli/fastpli/wiki/Solver): solve collisions of nerve fiber models
- [Simulation](https://github.com/3d-pli/fastpli/wiki/Simulation): 3D-PLI simulations of nerve fiber models

In addition, other modules exist to support io and analysis.

## Content

| module                  | information                                                   |
| ----------------------- | ------------------------------------------------------------- |
| `fastpli.analysis`      | analysis of 3D-PLI results                                    |
| `fastpli.io`            | input/output functions, e.g. to read/save fiber_bundles data  |
| `fastpli.model.sandbox` | building of simple 3d nerve fiber models                      |
| `fastpli.model.solver`  | generation of non intersection nerve fiber models             |
| `fastpli.objects`       | manipulation of fastpli objects (e.g. rotation)               |
| `fastpli.tools`         | mathematical tools and helper function                        |
| `fastpli.simulation`    | simulation of fiber models inside a virtual 3D-PLI microscope |

## Additional information

- [Wiki](https://github.com/3d-pli/fastpli/wiki)
- [Module documentation](https://3d-pli.github.io/fastpli/)
- [Examples and Jupyter notebooks](https://github.com/3d-pli/fastpli/tree/master/examples)

# Installation

Note: The current version of `fastpli` can only be executed under Linux as operating system because of dependencies.
If you want to use `fastpli` under Windows, please use the Windows Subsystem for Linux.
To allow for graphical output you have to install a X server.
Additional information can be found at <https://wiki.ubuntu.com/WSL>.

Support for macOS is planned for the future.

## Dependencies

### Requirements

- C++17
- Make
- CMake
- Python3
- MPI
- OpenGL (optional, recommended)

### Submodules

- pybind11

## Install instructions

### Packages

Install all necessary packages.

Example using Ubuntu or Debian:

```sh
sudo apt install gcc g++ cmake make git
sudo apt install python3-dev python3-venv
sudo apt install libopenmpi-dev freeglut3-dev
```

### Clone repository

```sh
git clone https://github.com/3d-pli/fastpli.git
cd fastpli
```

### Compilation

```sh
make fastpli
pip3 install .
```

# Examples

## Interactive jupyter notebooks

```sh
sandbox.ipynb
solver.ipynb
```

## Scripts

```sh
# install required modules for examples
pip3 install -r examples/requirements.txt

# run examples
python3 examples/sandbox.py
python3 examples/solver.py
python3 examples/simpli.py
python3 examples/simulation_pipeline.py

# run complete pipeline
python3 examples/crossing.py
```

# Tests

```sh
python3 setup.py test
```

# About this Project

## Libraries

All computationally intensive calculations are optimized either with **numba** on the Python side or with multithreading **c++**, which can be accessed via **pybind11**.
Additionally the simulation module supports the **Message Passing Interface (MPI)**.

## Contributions and Bug Reports

Please submit [issues](https://github.com/3d-pli/fastpli/issues) on GitHub to report
problems or suggest features. [Pull requests](https://github.com/3d-pli/fastpli/pulls)
are also welcome to add features or correct problems.
Please run the local env-CI environment `./CI/run-all.sh` or docker container `make docker` in advance.

## Literature

- [3D-PLI](https://dx.doi.org/10.3389%2Ffninf.2011.00034)
- [dense fiber modeling](https://arxiv.org/abs/1901.10284)
- [MEDUSA](https://doi.org/10.1016/j.neuroimage.2019.02.055)
- [simulation](https://doi.org/10.1016/j.neuroimage.2015.02.020)
- [tilting analysis](https://doi.org/10.3389/fnana.2018.00075)

## Authors

- **Felix Matuschke**

## References

|                                                                                                                                                                                  |                                                                                                                                                              |
| :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: | ------------------------------------------------------------------------------------------------------------------------------------------------------------ |
|                  [![Forschungszentrum Jülich](https://www.fz-juelich.de/SharedDocs/Bilder/INM/INM-1/EN/FZj_Logo.jpg?__blob=normal)](https://www.fz-juelich.de)                   | [Forschungszentrum Jülich](https://www.fz-juelich.de)                                                                                                        |
| [![FA-INM-1](https://avatars2.githubusercontent.com/u/51479655?s=200&v=4)](https://www.fz-juelich.de/inm/inm-1/EN/Forschung/Fibre%20Architecture/Fibre%20Architecture_node.html) | [Fiber Architecture - INM1 - Forschungszentrum Jülich](https://www.fz-juelich.de/inm/inm-1/EN/Forschung/Fibre%20Architecture/Fibre%20Architecture_node.html) |
|                                 [![HBP](https://sos-ch-dk-2.exo.io/public-website-production/img/HBP.png)](https://www.humanbrainproject.eu/en/)                                 | [Human Brain Project](https://www.humanbrainproject.eu/en/)                                                                                                  |

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
