
back to [readme](README.md)
# TODOs:
## main structure:
- [ ] analysis
  - [ ] epa
  - [ ] rofl
- [ ] model
  - [ ] generator
  - [x] solver
    - [ ] freeze objects
  - [ ] visualizer
    - [x] static
    - [ ] interactive
- [ ] objects
  - [x] fiber
  - [ ] tissue_container
- [ ] simulation
  - [ ] generator
  - [ ] simulator
  - [ ] apply_optic
- [ ] tools
  - [ ] h5io
  - [x] rotations
  
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
 - [ ] fiber
   - [ ] fiber::Data
   - [ ] fiber::Geometry
   - [ ] fiber::Segment (cone?)
   - [ ] fiber::Bundle
   - [ ] fiber::layer::Property
 - [ ] 