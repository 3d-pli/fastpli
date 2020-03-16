
[back](README.md)

# KNOWN ISSUES:
* multiprocessing has initial high cpu values with model.solver! 
* simpli.resolution does not call property
* nan values detected in light signal (for 90degree fibers?)
* simpli simulation has high initial copying and allocation time

# TODOs:
- [ ] consistent VOIs
- [ ] user friendly sandbox
- [ ] examples
- [ ] std::vector(float) for fibers
- [ ] helper_math.hpp
## simpli
- [ ] interpolation for mu, dn, ...
- [ ] step_size == dim.z
- [ ] rotate volume in simpli.generation
- [x] non untilt case
- [x] PM inverse light direction
- [ ] stokes vs jones
- [x] polarization value px, py
- [ ] check filter rotations
- [ ] check tilt direction
- [ ] check ALL input variables
- [ ] enlarge simulation for tilting <-> crop results
- [ ] np.array for cells
- [ ] optic: new_size = np.array(np.array(image.shape) // resize, dtype=int)
## VCS
- [ ] warning if fiber radius, segment length and r_min are questionable
- [ ] warning if collision close to 0 but "never ending"
- [ ] memory warning for fiber_bundles (gets quite big for small segment lengths)
  - [ ] non linear splitting and merging
- [ ] split volume first into equal cubes, and then into a octtree
- ideas:
- [ ] not moving Regions and moving Regions
- [ ] alowing only specific points of fibers to move
- [ ] alow movement anly along a specified axis
- [ ] move fibers as a group of bundles
- [ ] looping fibers
