# TODO

## Issues

* multiprocessing has initial high cpu values with model.solver
* simpli.pixel_size does not call property
* nan values detected in light signal (for 90degree fibers?)
* simpli simulation has high initial copying and allocation time

## simpli

* [ ] interpolation for mu, dn, ...
* [ ] stokes vs jones
* [ ] polarization value px, py
* [ ] test cells
* [ ] test mpi on multiple nodes

## VCS

* [ ] warning if fiber radius, segment length and r_min are questionable
* [ ] warning if collision close to 0 but "never ending"
* [ ] memory warning for fiber_bundles (gets quite big for small segment lengths)
  * [ ] non linear splitting and merging
* [ ] segment length per fiber, automatic determination for each fiber fun(f_radius)

---

[back](README.md)
