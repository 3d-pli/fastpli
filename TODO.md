
back to [readme](README.md)

# KNOWN ISSUES:
* memory usage of simpli unknown if voi is set before pixel_size
* simpli.resolution does not call property
* nan values detected in light signal (for 90degree fibers?)

# TODOs:
- [ ] std::vector(float) for fibers
- [ ] helper_math.hpp
## simpli
- [ ] check filter rotations
- [ ] check tilt direction
- [ ] check ALL input variables
- [ ] enlarge simulation for tilting <-> crop results
- [ ] too many arguments as property
- [ ] np.array for cells
- [ ] optic: new_size = np.array(np.array(image.shape) // resize, dtype=int)
- [ ] multisampling 
## VCS
- [ ] memory warning for fiber_bundles (gets quite big for small segment lengths)
  - [ ] non linear splitting and merging
- [ ] split volume first into equal cubes, and then into a octtree

## main structure:
- [x] analysis
  - [x] epa
  - [x] rofl
  - [ ] opt_pixel_resolution <- redo
- [ ] model
  - [ ] generator (FAconstructor)
  - [x] solver
    - [ ] freeze objects
  - [ ] visualizer
    - [x] static
    - [x] red for collision
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
