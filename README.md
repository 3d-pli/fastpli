<!-- 
________             ___________________________
___  __/_____ _________  /___  __ \__  /____  _/
__  /_ _  __ `/_  ___/  __/_  /_/ /_  /  __  /  
_  __/ / /_/ /_(__  )/ /_ _  ____/_  /____/ /   
/_/    \__,_/ /____/ \__/ /_/     /_____/___/    
-->
![](logo.png)

# Fiber Architecture Simulation Toolbox for PLI

## Installation:
```sh
# Compiling the source code and generating setup.py
make build

# installation with pip
pip3 install build/.
```

## examples:
```sh
# install required modules for examples
pip3 install -r examples/requirements.txt

# run examples
python3 examples/sandbox.py
python3 examples/model_solver.py
python3 examples/simpli.py
python3 examples/simulation_pipeline.py
```

# Dependencies
## Requirements:
 - C++
 - Make
 - CMake
 - Python3
 - MPI

## Optional Packages:
 - OpenGL
 - CUDA

## Submodules:
 - pybind11

# MPI execution:
The PLI simulation library simPLI supports mpi.

## Volume generation:
The generated volume will be splitt corresponding to the mpi processes. Each mpi process calculates the label_field and vector_field volume seperatly.

## PLI simulation:
Each mpi process simulates the light-tissue interaction on its volume. When a light beam needs to change to a neighboured volume, it is transmitted via mpi communication.

## Parallel HDF5io:
To be able to save the seperated data into a single hdf5 file, a paralel implementation of h5py is needed. However, a parallel installation of h5py is required. However, h5py can only be installed either serially or in parallel:
```sh
# serial:
make h5py-serial

# parallel:
make h5py-mpi

# clean:
h5py-clean
```

## example:
```sh
# simpli supports mpi 
mpiexec -n 2 python3 examples/simpli_mpi.py
```

# Additional Informations:
[program structure](docs/structure.md)

## TODOs:
[TODOs](docs/TODO.md)


---
## Authors
* **Felix Matuschke**: INM1 - Forschungszentrum Jülich


## License
This project is licensed under the GPLv3 License - see the [LICENSE](LICENSE) file for details
