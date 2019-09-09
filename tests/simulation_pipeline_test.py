import h5py
import numpy as np
import sys
import os

# fastpli
import fastpli.simulation
import fastpli.analysis

np.random.seed(42)

### This code is a simple reproducable test between version changes ###

if __name__ == '__main__':

    # PliGeneration ###
    simpli = fastpli.simulation.Simpli()
    simpli.pixel_size = 32
    simpli.resolution = 64
    simpli.set_voi([0, 640, 0, 640, 0, 64])
    fiber_radius = 64
    fiber_bundles = [[]]

    # corner
    fiber_bundles[-1].append(
        np.array([[0, 0, 30, 128],
                  [
                      simpli.dim[0] * simpli.pixel_size,
                      simpli.dim[1] * simpli.pixel_size, 30, fiber_radius
                  ]]))

    # left right up
    fiber_bundles[-1].append(
        np.array(
            [[0, simpli.dim[1] * simpli.pixel_size * 0.5, 30, fiber_radius],
             [
                 simpli.dim[0] * simpli.pixel_size,
                 simpli.dim[1] * simpli.pixel_size * 0.5, 30, fiber_radius
             ], [simpli.dim[0] * simpli.pixel_size, 0, 30, fiber_radius]]))

    # circle
    t = np.linspace(0, 2 * np.pi, 50)
    x = (np.sin(t) *
         (simpli.dim[0] // 3) + simpli.dim[0] // 2) * simpli.pixel_size
    y = (np.cos(t) *
         (simpli.dim[1] // 3) + simpli.dim[1] // 2) * simpli.pixel_size
    z = t * 0 + simpli.dim[2] * simpli.pixel_size / 2
    r = t * 0 + fiber_radius
    fiber_bundles[-1].append(np.array([x, y, z, r]).T)

    simpli.fiber_bundles = fiber_bundles
    simpli.fiber_bundles_properties = [[(0.5, 0.004, 1, 'p'),
                                        (1.0, 0.004, 1, 'r')]]

    simpli.memory_usage()
    label_field, vec_field, tissue_properties = simpli.generate_tissue()

    # PliSimulation ###
    simpli.filter_rotations = np.deg2rad([90, 120, 150, 0, 30, 60])
    simpli.light_intensity = 26000
    simpli.untilt_sensor = True
    simpli.wavelength = 525

    tilts = [(0, 0), (5.5, 0), (5.5, 90), (5.5, 180), (5.5, 270)]
    image_stack = []

    # print("Run Simulations:")
    for t, (theta, phi) in enumerate(tilts):
        print('Step:', t, theta, phi)
        image = simpli.run_simulation(label_field, vec_field, tissue_properties,
                                      np.deg2rad(theta), np.deg2rad(phi))
        image = simpli.apply_optic(image)

        epa = simpli.apply_epa(image)

        image_stack.append(image)

    # save mask for analysis, but restore None value so that everything gets analysed
    mask = np.sum(label_field, 2) > 0
    mask = simpli.apply_resize_mask(mask) > 0.1

    # print("Run ROFL:")
    rofl_direction, rofl_incl, rofl_t_rel, _ = simpli.apply_rofl(image_stack,
                                                                 mask=mask)

    # print("Unit Vectors")
    unit_x, unit_y, unit_z = fastpli.analysis.images.unit_vectors(
        rofl_direction, rofl_incl, mask)

    # print("FOMs")
    img = fastpli.analysis.images.fom_hsv_black(rofl_direction, rofl_incl, mask)
    img = fastpli.analysis.images.fom_rgb(rofl_direction, rofl_incl, mask)
    img = fastpli.analysis.images.hsvblack_sphere(n=64)
    img = fastpli.analysis.images.rgb_sphere(n=64)
