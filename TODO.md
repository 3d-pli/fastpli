
back to [readme](README.md)
# TODOs:
- [ ] git submodule to cmake
- [ ] 'make build' output to fastpli
- [ ] simpli.voxel_size
- [ ] check filter rotations
- [ ] check tilt direction
- [ ] enlarge simulation for tilting <-> crop results

## main structure:
- [x] analysis
  - [x] epa
  - [x] rofl
  - [ ] opt_pixel_resolution
- [ ] model
  - [ ] generator (FAconstructor)
  - [x] solver
    - [ ] freeze objects
  - [ ] visualizer
    - [x] static
    - [ ] red for collision
    - [ ] interactive
- [x] objects
  - [x] fiber
  - [x] tissue_container
  - [x] tissue_container as contigious np.array without copy
- [ ] simulation
  - [x] generator
  - [x] simulator
    - [ ] PM light direction 
    - [ ] rotation angles
    - [ ] stokes vs jones
    - [ ] polarisation value px, py
  - [x] *.txt output for label_field
  - [x] apply_optic
- [ ] tools
  - [ ] h5io
    - [ ] mpi - 4gb threshold
  - [x] rotation

## examples:
- [ ] model
- [ ] simpli
- [ ] pipeline (generation, simulation, analysis)
  
## comparison to older programs:
- [ ] simulation <-> simPLI
- [ ] model.solver <-> vcs

## HPC:
- [ ] JURON
- [x] JURECA

## Tests
- [x] CI
- [x] docker-file
- [ ] simulation.generator
- [ ] simulation.simulator

## restructure
 - [x] fiber
   - [x] fiber::Fiber
   - [x] fiber::Geometry (model)
   - [x] fiber::layer::Property (simulation)
   - [x] fiber::Bundle (simulation)

## include:
- [x] aabb.hpp
- [x] vemath.hpp

# File Structure:
## PyPackage structure:
```sh
`-- fastpli
    |-- analysis
    |-- model
    |   |-- generator
    |   |-- solver
    |   `-- visualizer
    |-- objects
    |-- simulation
    `-- tools
        `-- h5io
```

## C++ structure:
### module represented in cpp:
```sh
|-- module-name
|   |-- CMakeLists.txt
|   |-- bindings
|   |   `-- foo_module.cpp
|   |-- foo.cpp
|   `-- foo.hpp

```

## Tests:
```sh
|-- tests
    `-- module-name
        `-- foo_test.py
```
