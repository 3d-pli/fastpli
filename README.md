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
make install
```

## Known Issues:
If you run a python script with 'mpirun' and use h5py, h5py has to be installed with 'make h5py-mpi'. Only one version of h5py can be installed.

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

## TODOs:
look [here](TODO.md)

## Authors
* **Felix Matuschke**

## License
This project is licensed under the GPLv3 License - see the [LICENSE](LICENSE) file for details
