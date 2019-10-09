import unittest

import numpy as np

import fastpli.simulation
import fastpli.analysis
import fastpli.tools
import fastpli.objects

from fastpli.analysis._ROFL_with_jacobi import _calc_Int_single_fiber_fitted


class MainTest(unittest.TestCase):

    def test_simple_rofl(self):

        simpli = fastpli.simulation.Simpli()
        simpli.omp_num_threads = 1
        simpli.voxel_size = 60  # in mu meter
        simpli.dim = [1, 1, 1]

        # single voxel
        label_field = np.ones((1), dtype=np.int32)
        vec_field = np.zeros((1, 3), dtype=np.float32)
        vec_field[:] = [np.cos(np.deg2rad(45)), 0, np.sin(np.deg2rad(45))]
        tissue_properties = np.array([[0, 0], [-0.001, 0]])

        # Simulate PLI Measurement ###
        simpli.filter_rotations = np.deg2rad([0, 30, 60, 90, 120, 150])
        simpli.light_intensity = 26000  # a.u.
        simpli.interpolate = True
        simpli.untilt_sensor_view = False
        simpli.wavelength = 525  # in nm
        simpli.step_size = 1  # in voxel_size
        TILTS = [(0, 0), (5.5, 0), (5.5, 90), (5.5, 180), (5.5, 270)]

        tilting_stack = [None] * len(TILTS)

        for t, (theta, phi) in enumerate(TILTS):
            simpli.step_size = 1 / np.cos(np.deg2rad(theta))
            images = simpli.run_simulation(label_field, vec_field,
                                           tissue_properties, np.deg2rad(theta),
                                           np.deg2rad(phi))

            # calculate modalities
            tilting_stack[t] = images

        rofl_direction, rofl_incl, rofl_t_rel, _ = simpli.apply_rofl(
            tilting_stack,
            tilt_angle=np.deg2rad(TILTS[-1][0]),
            gain=3,  # only LAP!
            dir_offset=0,
            num_threads=1,
            grad_mode=False)

        t_rel = 4 * simpli.voxel_size * abs(
            tissue_properties[1][0]) / (simpli.wavelength / 1e3)

        self.assertTrue(
            abs(rofl_direction.flatten()[0] - np.deg2rad(0)) < 1e-10 or
            abs(rofl_direction.flatten()[0] - np.deg2rad(180)) < 1e-10)
        self.assertTrue(
            abs(rofl_incl.flatten()[0] - np.deg2rad(45) < 1e-9) or
            abs(rofl_incl.flatten()[0] - np.deg2rad(-45)) < 1e-9)
        self.assertTrue(abs(rofl_t_rel.flatten()[0] - t_rel) < 1e-8)


if __name__ == '__main__':
    unittest.main()
