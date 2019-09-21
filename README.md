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
### with venv:
```sh
make install
```
### with pip:
```sh
make build
pip3 install build/
```

### examples:
```sh
pip3 install -r examples/requirements.txt
```

## Parallel HDF5io:
h5py has to be installed with mpi flags: 
```sh
make h5py-mpi
```
Only one version of h5py can be installed.

## Required Programs:
 - G++/Clang
 - Make
 - CMake
 - Python3
 - HDF5
 - MPI
 - OpenMP

## Optional Packages:
 - OpenGL
 - CUDA

## Used Libraries:
 - pybind11

## Additional Informations:
[program structure](docs/structure.md)

## TODOs:
[TODOs](docs/TODO.md)

## Authors
* **Felix Matuschke**

## License
This project is licensed under the GPLv3 License - see the [LICENSE](LICENSE) file for details
