
back to [readme](README.md)
# TODOs:
## main structure:
- [ ] analysis
  - [ ] epa
  - [ ] rofl
- [ ] model
  - [ ] generator (FAconstructor)
  - [x] solver
    - [ ] freeze objects
  - [ ] visualizer
    - [x] static
    - [ ] interactive
- [ ] objects
  - [x] fiber
  - [x] tissue_container
  - [ ] tissue_container as contigious np.array without copy
- [ ] simulation
  - [x] generator
  - [x] simulator
  - [ ] apply_optic
- [ ] tools
  - [ ] h5io
  - [x] rotation
  
## comparison to older programs:
- [ ] simulation <-> simPLI
- [ ] model.solver <-> vcs

## HPC:
- [ ] JURON
- [ ] JURECA

## Tests
- [x] CI
- [x] docker-file

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

## include:
 - aabb.hpp
 - vemath.hpp

## restructure
 - [x] fiber
   - [x] fiber::Fiber
   - [x] fiber::Geometry (model)
   - [x] fiber::layer::Property (simulation)
   - [x] fiber::Bundle (simulation)

