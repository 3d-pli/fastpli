<!-- 
________             ___________________________
___  __/_____ _________  /___  __ \__  /____  _/
__  /_ _  __ `/_  ___/  __/_  /_/ /_  /  __  /  
_  __/ / /_/ /_(__  )/ /_ _  ____/_  /____/ /   
/_/    \__,_/ /____/ \__/ /_/     /_____/___/    
-->
![](logo.png)

# Fiber Architecture Simulation Toolbox for PLI

## Basic Information:
`fastpli` is a Python package consisting of the following modules

| module                  | information                                                  |
| ----------------------- | ------------------------------------------------------------ |
| `fastpli.analysis`      | analysis of 3D-PLI results                                   |
| `fastpli.io`            | input/output functions, e.g. to read/save fiber_bundles data |
| `fastpli.model.sandbox` | building of simple 3d nerve fiber models                     |
| `fastpli.model.solver`  | generation of non intersection nerve fiber models            |
| `fastpli.objects`       | manipulation of fastpli objects (e.g. rotation)              |
| `fastpli.tools`         | mathematical tools and helper function                       |
| `fastpli.simulation`    | simulation of fiber models inside a virtual 3D-PLI microscop |

The aim of this package is to provide consistent system that allows the following: 
* **model** 3d (non colliding) nerve fibers
* **simulate** nerve fiber inside a virtual tiltable 3D-PLI microscop
* **analyse** the simulated signals to extract the resulting fiber orientation

See **Wiki** for detailed information

## Performance
All computationally intensive calculations are optimized either with **numba** on the Python side or with multithreading **c++**, which can be accessed via **pybind11**. Additionally the simulation module supports the **Message Passing Interface (MPI)**.

# Installation:

## Dependencies
### Requirements:
 - C++17
 - Make
 - CMake
 - Python3
 - MPI
 - OpenGL

### Submodules:
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
<!-- libhdf5-openmpi-dev -->

Example using Archlinux:

```
pacman --noconfirm --needed -Syu gcc make cmake git
pacman --noconfirm --needed -Syu python python-pipenv python-pip
pacman --noconfirm --needed -Syu openmpi
pacman --noconfirm --needed -Syu freeglut glu
```
<!--  hdf5-openmpi -->

### Clone repository

```sh
git clone git@github.com:3d-pli/fastPLI.git
cd fastpli
```

### Compilation 

#### Makefile

The Makefile contains the instruction to create a virtual Python environment `env`. 
All the following instructions and examples use this virtual environment.

```sh
make install
```

#### Or Step by Step

The Makefile runs the following commands (with a creation of an virtual python environment):

```sh
git submodule update --init
mkdir build
cp -al --remove-destination src/fastpli build/
cd build/
cmake ..
make -j
pip3 install .
```

It creates a hardlink so that a pips debugging process is available.
For the future it is planned that the compilation process will be handled by pip.
However, the standard approach of this process has in the past significantly slowed down the compilation time.

## Running examples:

```sh
# install required modules for examples
env/bin/pip3 install -r examples/requirements.txt

# run examples
env/bin/python3 examples/sandbox.py
env/bin/python3 examples/solver.py
env/bin/python3 examples/simpli.py
env/bin/python3 examples/simulation_pipeline.py
```

## Authors
* **Felix Matuschke**

## References
|                                                                                                                                                                                                                |                                                                                                                                                              |
| :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: | ------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| [![](https://www.fz-juelich.de/SharedDocs/Bilder/INM/INM-1/DE/PLI/PLI-GruppenLogo.png?__blob=thumbnail)](https://www.fz-juelich.de/inm/inm-1/EN/Forschung/Fibre%20Architecture/Fibre%20Architecture_node.html) | [Fiber Architecture - INM1 - Forschungszentrum Jülich](https://www.fz-juelich.de/inm/inm-1/EN/Forschung/Fibre%20Architecture/Fibre%20Architecture_node.html) |
|                                                 [![](https://sos-ch-dk-2.exo.io/public-website-production/img/HBP.png)](https://www.humanbrainproject.eu/en/)                                                  | [Human Brain Project](https://www.humanbrainproject.eu/en/)                                                                                                  |

## License
This project is licensed under the GPLv3 License - see the [LICENSE](LICENSE) file for details
