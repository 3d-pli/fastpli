<!--
________             ___________________________
___  __/_____ _________  /___  __ \__  /____  _/
__  /_ _  __ `/_  ___/  __/_  /_/ /_  /  __  /
_  __/ / /_/ /_(__  )/ /_ _  ____/_  /____/ /
/_/    \__,_/ /____/ \__/ /_/     /_____/___/
-->

# Fiber Architecture Simulation Toolbox for 3D-PLI

![fastpli-logo](logo.svg)

The
[Fiber Architecture Simulation Toolbox for 3D-PLI (fastpli)](https://github.com/3d-pli/fastpli)
is a toolbox for
[polarized light imaging (PLI)](https://www.fz-juelich.de/inm/inm-1/EN/Forschung/Fibre%20Architecture/Fibre%20Architecture_node.html)
with three main purposes:

<img align="right" src="https://raw.githubusercontent.com/wiki/3d-pli/fastpli/images/fiber.png" alt="fiber" width="125">

- [`Sandbox` - designing of nerve fiber models](https://github.com/3d-pli/fastpli/wiki/NerveFiber):
  The first module allows the user to create different types of nerve fiber
  bundles and additionally fill them with individual nerve fibers.

  - [Details](https://github.com/3d-pli/fastpli/wiki/NerveFiber)
  - [Tutorial](https://github.com/3d-pli/fastpli/wiki/tutorial-sandbox)

<img align="right" src="https://raw.githubusercontent.com/wiki/3d-pli/fastpli/images/solver_1_cropped.gif" alt="fiber" width="125">

- [`Solver` - generating collision free models](https://github.com/3d-pli/fastpli/wiki/Solver):
  The second module takes as input a configuration of nerve fibers and checks
  them for spatial collisions. Since nerve fibers cannot overlap in reality, one
  must ensure that the models follow the same rules. The solver module
  implements a simple algorithm that checks for collisions and, if it finds any,
  pushes the colliding segments of the fibers slightly apart. This is repeated
  until all collisions are solved.

  - [Details](https://github.com/3d-pli/fastpli/wiki/Solver)
  - [Tutorial](https://github.com/3d-pli/fastpli/wiki/tutorial-solver)

<img align="right" src="https://raw.githubusercontent.com/wiki/3d-pli/fastpli/images/optic_chiasm_10_0.png" alt="fiber" width="125">

- [`Simulation` - simulation of 3D-Polarized Light Imaging](https://github.com/3d-pli/fastpli/wiki/Simulation):
  The simulation module enables the simulation of
  [3D Polarized Light Imaging (3D-PLI)](https://www.fz-juelich.de/inm/inm-1/EN/Forschung/Fibre%20Architecture/Fibre%20Architecture_node.html).
  This is a microscopic technique that allows the polarization change of light
  moving through a brain section to be measured. Due to the birefringence
  property of the myelin surrounding the nerve fibers, the polarization state
  changes. This change enables the calculation of the 3d orientation of the
  nerve fibers in the brain slice.

  - [Details](https://github.com/3d-pli/fastpli/wiki/Simulation)
  - [Tutorial](https://github.com/3d-pli/fastpli/wiki/tutorial-simulation)

## Wiki

<https://github.com/3d-pli/fastpli/wiki>

## Example

As an example, a simplified model of the optic chiasm is presented. This
structure in the brain allows nerve fibers from the eyes to cross each other and
connect to the opposite side of the brain. In addition, a certain portion
remains on the same side of the brain.

- [Tutorial](https://github.com/3d-pli/fastpli/wiki/tutorial-optic_chiasm)

![png](https://raw.githubusercontent.com/wiki/3d-pli/fastpli/optic_chiasm_files/optic_chiasm_18_1.png)
![png](https://raw.githubusercontent.com/wiki/3d-pli/fastpli/optic_chiasm_files/optic_chiasm_19_0.png)

## Module lists

[API documentation](https://3d-pli.github.io/fastpli/)

| module                                                                                               | information                                                   |
| ---------------------------------------------------------------------------------------------------- | ------------------------------------------------------------- |
| [`fastpli.analysis` ](https://3d-pli.github.io/fastpli/_autosummary/fastpli.analysis.html)           | analysis of 3D-PLI results                                    |
| [`fastpli.io` ](https://3d-pli.github.io/fastpli/_autosummary/fastpli.io.html)                       | input/output functions, e.g. to read/save fiber_bundles data  |
| [`fastpli.model.sandbox` ](https://3d-pli.github.io/fastpli/_autosummary/fastpli.model.sandbox.html) | building of simple 3d nerve fiber models                      |
| [`fastpli.model.solver` ](https://3d-pli.github.io/fastpli/_autosummary/fastpli.model.solver.html)   | generation of non intersection nerve fiber models             |
| [`fastpli.objects` ](https://3d-pli.github.io/fastpli/_autosummary/fastpli.objects.html)             | manipulation of fastpli objects (e.g. rotation)               |
| [`fastpli.tools` ](https://3d-pli.github.io/fastpli/_autosummary/fastpli.tools.html)                 | mathematical tools and helper function                        |
| [`fastpli.simulation` ](https://3d-pli.github.io/fastpli/_autosummary/fastpli.simulation.html)       | simulation of fiber models inside a virtual 3D-PLI microscope |

# Installation

### Note:

> The current version of `fastpli` can only be run under Linux as operating
> system due to dependencies. If you want to use `fastpli` under Windows, please
> use the Windows subsystem for Linux. To enable graphical output, you must
> install an X server. For more information, see <https://wiki.ubuntu.com/WSL>.
> Support for macOS is planned for the future.

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

For Ubuntu:

```sh
sudo apt update
sudo apt install gcc g++ cmake make git
sudo apt install python3-dev python3-venv
sudo apt install libopenmpi-dev freeglut3-dev
```

### Clone repository

```sh
git clone --recursive https://github.com/3d-pli/fastpli.git
cd fastpli
```

### Compilation

Use your favorite environment e. g. `python3 -m venv env` and
`source env/bin/activate`. Update your pip version with `pip3 install pip -U`.

```sh
make fastpli
pip3 install .
```

# Examples

## Tutorials

```sh
# install required modules for examples
pip3 install -r examples/requirements.txt

jupyter-notebook examples/sandbox.ipynb
jupyter-notebook examples/solver.ipynb
jupyter-notebook examples/simulation.ipynb
jupyter-notebook examples/optic_chiasm.ipynb
```

## Scripts

```sh
# install required modules for examples
pip3 install -r examples/requirements.txt

# run examples
python3 examples/sandbox.py
python3 examples/solver.py
python3 examples/simulation.py
python3 examples/optic_chiasm.py
```

# Tests

```sh
python3 setup.py test
```

# About this Project

## Libraries

All computationally intensive calculations are optimized either with **numba**
on the Python side or with multithreading **C++**, which can be accessed via
**pybind11**. Additionally the simulation module supports the **Message Passing
Interface (MPI)**.

## Contributions and Bug Reports

Please submit [issues](https://github.com/3d-pli/fastpli/issues) on GitHub to
report problems or suggest features.
[Pull requests](https://github.com/3d-pli/fastpli/pulls) are also welcome to add
features or correct problems. Please run the local env-CI environment
`./CI/run-all.sh` or docker container `make docker` in advance.

## Literature

- [3D-PLI](https://dx.doi.org/10.3389%2Ffninf.2011.00034)
- [dense fiber modeling](https://arxiv.org/abs/1901.10284)
- [MEDUSA](https://doi.org/10.1016/j.neuroimage.2019.02.055)
- [simulation](https://doi.org/10.1016/j.neuroimage.2015.02.020)
- [tilting analysis](https://doi.org/10.3389/fnana.2018.00075)

## Authors

- **Felix Matuschke**

## References

[fastPLI](https://github.com/3d-pli/fastpli) is an open source toolbox for
modeling nerve fibers, simulating them in a
[3D-PLI](https://dx.doi.org/10.3389%2Ffninf.2011.00034) microscope and the
signal processing developed by the
[fiber architecture group](https://www.fz-juelich.de/inm/inm-1/EN/Forschung/Fibre%20Architecture/Fibre%20Architecture_node.html)
at the [Forschungszentrum Jülich](https://www.fz-juelich.de) -
[INM1](https://www.fz-juelich.de/inm/inm-1/EN/Home/home_node.html). This project
has received funding from the European Union’s Horizon 2020 Research and
Innovation Programme under Grant Agreement No. 7202070
([Human Brain Project](https://www.humanbrainproject.eu/en/) SGA2).

|                                                                                                                                                                                  |                                                                                                                                                              |
| :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: | ------------------------------------------------------------------------------------------------------------------------------------------------------------ |
|       [![Forschungszentrum Jülich](https://upload.wikimedia.org/wikipedia/commons/4/40/Logo_des_Forschungszentrums_J%C3%BClich_seit_2018.svg)](https://www.fz-juelich.de)        | [Forschungszentrum Jülich](https://www.fz-juelich.de)                                                                                                        |
| [![FA-INM-1](https://avatars2.githubusercontent.com/u/51479655?s=200&v=4)](https://www.fz-juelich.de/inm/inm-1/EN/Forschung/Fibre%20Architecture/Fibre%20Architecture_node.html) | [Fiber Architecture - INM1 - Forschungszentrum Jülich](https://www.fz-juelich.de/inm/inm-1/EN/Forschung/Fibre%20Architecture/Fibre%20Architecture_node.html) |
|                                 [![HBP](https://sos-ch-dk-2.exo.io/public-website-production/img/HBP.png)](https://www.humanbrainproject.eu/en/)                                 | [Human Brain Project](https://www.humanbrainproject.eu/en/)                                                                                                  |

## License

This project is licensed under the MIT License - see the
[LICENSE](https://github.com/3d-pli/fastpli/blob/main/LICENSE) file for details
