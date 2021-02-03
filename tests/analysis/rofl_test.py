import unittest

import numpy as np

import fastpli.simulation
import fastpli.analysis
import fastpli.tools
import fastpli.objects

# TODO: test rofl, not simpli.apply_rofl


class MainTest(unittest.TestCase):
    def test_simple_rofl(self):

        simpli = fastpli.simulation.Simpli()
        simpli.omp_num_threads = 1
        simpli.voxel_size = 60  # in micro meter
        simpli.dim = [1, 1, 1]

        # single voxel
        tissue = np.ones((1, 1, 1), dtype=np.int32)
        optical_axis = np.zeros((1, 1, 1, 3), dtype=np.float32)
        optical_axis[:] = [np.cos(np.deg2rad(45)), 0, np.sin(np.deg2rad(45))]
        tissue_properties = np.array([[0, 0], [-0.001, 0]])

        # Simulate PLI Measurement ###
        simpli.filter_rotations = np.deg2rad([0, 30, 60, 90, 120, 150])
        simpli.light_intensity = 26000  # a.u.
        simpli.interpolate = "Slerp"
        simpli.untilt_sensor_view = False
        simpli.wavelength = 525  # in nm
        simpli.step_size = 1  # in voxel_size
        simpli.tilts = np.deg2rad([(0, 0), (5.5, 0), (5.5, 90), (5.5, 180),
                                   (5.5, 270)])
        tilting_stack = [None] * 5

        for t, (theta, phi) in enumerate(simpli.tilts):
            simpli.step_size = 1 / np.cos(theta)
            images = simpli.run_simulation(tissue, optical_axis,
                                           tissue_properties, theta, phi)

            # calculate modalities
            tilting_stack[t] = images

        simpli.noise_model = lambda x: np.random.negative_binomial(
            x / (3 - 1), 1 / 3)
        rofl_direction, rofl_incl, rofl_t_rel, _ = simpli.apply_rofl(
            tilting_stack, grad_mode=False)

        t_rel = 4 * simpli.voxel_size * abs(
            tissue_properties[1][0]) / (simpli.wavelength / 1e3)

        self.assertTrue(
            abs(rofl_direction.flatten()[0] - np.deg2rad(0)) < 1e-10
            or abs(rofl_direction.flatten()[0] - np.deg2rad(180)) < 1e-10)
        self.assertTrue(
            abs(rofl_incl.flatten()[0] - np.deg2rad(45) < 1e-9)
            or abs(rofl_incl.flatten()[0] - np.deg2rad(-45)) < 1e-9)
        self.assertTrue(abs(rofl_t_rel.flatten()[0] - t_rel) < 1e-8)


if __name__ == '__main__':
    unittest.main()
