
[readme](README.md)

# KNOWN ISSUES:
* multiprocessing has initial high cpu values with model.solver! 
* simpli.resolution does not call property
* nan values detected in light signal (for 90degree fibers?)

# TODOs:
- [ ] consistent vois
- [ ] user friendly sandbox
- [ ] examples
- [ ] std::vector(float) for fibers
- [ ] helper_math.hpp
## simpli
- [ ] interpolation for mu, dn, ...
- [ ] stepsize == dim.z
- [ ] rotate volume in simipli.generation
- [x] non untilt case
- [x] PM inverse light direction
- [ ] stokes vs jones
- [x] polarisation value px, py
- [ ] check filter rotations
- [ ] check tilt direction
- [ ] check ALL input variables
- [ ] enlarge simulation for tilting <-> crop results
- [ ] np.array for cells
- [ ] optic: new_size = np.array(np.array(image.shape) // resize, dtype=int)
## VCS
- [ ] warning if fiber radius, segmentlength and rmin are questionable
- [ ] warning if collision close to 0 but "never ending"
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
