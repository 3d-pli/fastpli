---
title: "fastPLI: A Fiber Architecture Simulation Toolbox for 3D-PLI"
tags:
  - 3D-PLI
  - Python
  - microscopy
  - simulation
authors:
  - name: Felix Matuschke
    orcid: 0000-0001-6435-5351
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Katrin Amunts
    orcid: 0000-0001-5828-0867
    affiliation: 1, 2
  - name: Markus Axer
    orcid: 0000-0001-5871-9331
    affiliation: 1
affiliations:
  - name: Institute of Neuroscience and Medicine (INM-1), Forschungszentrum Jülich GmbH, 52425, Jülich, Germany
    index: 1
  - name: Cécile and Oskar Vogt Institute for Brain Research, University Hospital Düsseldorf, University of Düsseldorf, 40204, Düsseldorf, Germany.
    index: 2
date: 18. December 2020
bibliography: paper.bib
---

# Statement of need

3D Polarized Light Imaging (3D-PLI) is a microscopic neuroimaging technique used to study the nerve fiber architecture in unstained histological brain sections at the micrometer scale [@Axer2011].
It provides image contrast for fibers and fiber tracts, and ultimately enables reconstruction of 3D nerve fiber orientations.
The physical effect behind 3D-PLI is the optical property of the nerve fibers called birefringence.
Due to this intrinsic birefringence, it is possible to use polarized light, pass it through a thin brain section and observe the change of the polarization state of light.
This change is directly related to the 3D orientation of the fibers and also provides strong contrasts between fibers and other tissue components.

To understand the influence of the underlying fiber structure, simulations are an essential tool. They allow testing different hypotheses, knowing the ground truth. It has already been shown that simulations with scattered light within tissue sections require models with irregularities to mimic the behavior of the scattered light. This knowledge can now be used to understand structures such as fiber crossings that have been difficult to interpret in 3D-PLI [@Menzel2020].

In addition, the generated nerve fiber models can be used in other imaging simulation techniques such as diffusion magnetic resonance imaging (dMRI).

In recent years, various software tools have been developed to design fibre models.
Nerve fiber modeling is commonly used in dMRI.
But many modeling techniques do not use volumetric representations, especially collision-free ones.
In the last decade, an increasing number of algorithms for non-trivial overlapping structures have been developed [@Altendorf2011; @Chapelle2015; @Mingasson2017; @Ginsburger2019].
While these algorithms are specialized in their field, _fastPLI_ also provides a dedicated tool for 3D-PLI simulation based on linear optics.
Here, the focus is on the previously developed algorithm [@Matuschke2019], which provides a fast method to generate collision-free results for white matter structures in the brain.

Different types of simulations for polarized light are for example described in [RamellaRoman2005; @vanTurnhout2009; @Jiang2020].
However, to our knowledge, none of these techniques have been used to simulate the effects of polarized light on nerve fibers, except `simPLI` [@Dohmen2015], which is included in this toolbox.

# Summary

_fastPLI_ is an open source toolbox based on Python and C++ for modeling myelinated axons, i.e. nerve fibers and simulating the results of measurement of fiber orientations with a polarization microscope using 3D-PLI.

The _fastPLI_ package includes the following modules:

1. **Fiber Modelling Modules:**
   A detailed 3D modelling of nerve fibers at the micrometer level is essential as input for the measurement simulation.
   In order to recreate biological tissue as a model, it is important that the nerve fibers do not spatially overlap.
   We have decided to implement a solver module that takes any configuration of fiber bundles as input and converts it over several iterations into a collision-free configuration.
   In order to generate collision free fiber arrangements, a dedicated algorithm to prohibit such overlaps has been developed [@Matuschke2019].

2. **Simulation Module:**
   The 3D-PLI simulation is based on Stokes vector and Müller matrix approaches as described in [@Dohmen2015; @Menzel2015].
   For the simulation the polarimetric setup can be equipped with a tiltable specimen stage [@Axer2011; @Schmitz2018].
   By this means the brain section can be scanned from oblique views which adds important information to unambiguously analyze the 3D fiber orientation.

3. **Analysis Module:**
   The resulting simulated measurements (i.e., image stacks of a section acquired at different polarizing filter rotation angles and, optionally, at different oblique views) can be processed similarly to the real, experimental 3D-PLI [@Axer2011; @Schmitz2018].

All computationally intensive calculations are optimized either with _numba_ on the Python side or with multithreading _C++_ algorithms, which can be accessed via _pybind11_ inside the Python package [@Lam2015;@pybind11].
Additionally, the simulation module supports the Message Passing Interface (MPI) to facilitate the simulation of very large volumes on multiple computer nodes.

# Installation

The _fastPLI_ package has to be built with a _C++17_ compiler.
A Makefile allows a simple local installation.
It generates the necessary libraries inside a build folder and a matching setup.py file, which can be used for a second installation process with pip.

```sh
make fastpli
pip3 install .
```

All necessary configurations are handled in the background with CMake.
All required software libraries are listed in the [software repository](https://github.com/3d-pli/fastpli).

# Usage & Examples

A more detailed description of _fastPLI_ including examples, jupyter notebooks and tutorials can be found in the software package [^1], at the wiki pages [^2] and the API documentation [^3].

[^1]: [https://github.com/3d-pli/fastpli](https://github.com/3d-pli/fastpli/tree/main/examples)
[^2]: [https://github.com/3d-pli/fastpli/wiki](https://github.com/3d-pli/fastpli/wiki)
[^3]: [https://3d-pli.github.io/fastpli/](https://3d-pli.github.io/fastpli/)

```sh
python3 examples/sandbox.py
python3 examples/solver.py
python3 examples/simulation.py
python3 examples/optic_chiasm.py
```

## Fiber Modelling

Two modules exist to allow the user to build non-colliding white matter nerve fiber models: `fastpli.model.sandbox` and `fastpli.model.solver`.

The `fastpli.model.sandbox` module contains multiple functions to build from coordinates single individual fibers up to fiber bundles composed of individual fibers.

The `fastpli.model.solver` module contains a `Solver` class to generate non-colliding fibers constellations from an initial input dataset [@Matuschke2019].

An example of the solving process is shown in fig. \ref{generation} a).
The red colored segments indicate that a collision with this segment is detected.
With further iterations, the number of collisions decreases.
At the end, a collision-free configurations results.
Figure \ref{generation} b) shows the cross-section through the resulting fiber configuration.
Each fiber bundle has a different gray value.

![a) Solving process of two crossing fiber bundles. Individual colliding segments are colored in red. b) Cross-section of the collision solving process.\label{generation}](generation.png){width=100%}

The resulting orientations for each segment can be then visualized in a polar histogram (see fig. \ref{orientation}).

![orientation distribution a) initial, b) resulting fiber configuration. \label{orientation}](orientation_plot.png){width=100%}

## Simulation and Analysis

The module `fastpli.simulation` contains the class `Simpli`.
This class contains all tools required to simulate the 3D-PLI setup [@Dohmen2015; @Menzel2015] and to analyze the generated images.
Final results include the stack of images as well as the derived modalities referred to as transmittance, direction, retardation, inclination, relative thickness, and 3D fiber orientation (FOM) maps [@Schmitz2018].

An example of the simulation results is shown in fig. \ref{crossingmodalities}.
a) indicates the resulting transmittance map.
It represents the mean light intensity transmission through the cut tissue.
b) represents the fiber direction in the tissue plane.
c) shows the resulting retardation, which is a measure of the amplitude of the resulting sinusoidal signal from the simulation in each pixel.

![PLI modalities: a) transmittance, b) direction, c) retardation\label{crossingmodalities}](crossing_modalities.png){width=100%}

The data can be further analyzed with an advanced tilting simulation and analysis \ref{crossingrofl}.
a) represents the fiber direction in the tissue plane.
b) visualizes the fiber inclination in the tissue plane.
c) shows the relative thickness of the fibers within a pixel volume.
d) represents the resulting FOM. The colors encode the 3D orientation.

![Tilting analysis results: a) direction, b) inclination, c) relative thickness, d) FOM\label{crossingrofl}](crossing_rofl.png){width=100%}

The analysis is accessible inside the simulation module via a _pipeline_ inside the `fastpli.simulation.Simpli` class. For further information see `examples/simulation_pipeline.py`.

# Acknowledgements

This study was funded by the European Union's Horizon 2020 Research and Innovation Programme under Grant Agreement No. 785907 (Human Brain Project SGA2) and No. 945539 (Human Brain Project SGA3).
We gratefully acknowledge the computing time granted through JARA-HPC on the supercomputer JURECA [@jureca] at Forschungszentrum Jülich, Germany.

# Author Contributions

F.M. developed the software and wrote the documentation.
M.A. provided guidance throughout the project.
K.A. and M.A. supervised the project.
F.M. wrote the manuscript with input from all authors.

# References
