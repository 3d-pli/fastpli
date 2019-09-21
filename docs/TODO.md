
back to [readme](README.md)

# KNOWN ISSUES:
* memory usage of simpli unknown if voi is set before pixel_size
* simpli.resolution does not call property
* nan values detected in light signal (for 90degree fibers?)

# TODOs:
- [ ] std::vector(float) for fibers
- [ ] helper_math.hpp
## simpli
- [ ] PM inverse light direction
- [ ] stokes vs jones
- [ ] polarisation value px, py
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
  
## HPC:
- [ ] JURON
- [x] JURECA

## Tests
- [x] CI
- [x] docker-file
- [ ] simulation.generator
- [ ] simulation.simulator
