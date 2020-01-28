import unittest
import numpy as np
import h5py
import os
import warnings

from fastpli.simulation import Simpli


class MainTest(unittest.TestCase):

    def setUp(self):
        self.fiber_bundles = [[[[0, 0, 0, 1], [1, 1, 1, 1], [2, 2, 2, 1]]]]
        self.fiber_bundles_properties = [[(0.333, -0.004, 10, 'p'),
                                          (0.666, -0.004, 5, 'b'),
                                          (1.0, 0.004, 1, 'r')]]

        self.simpli = Simpli()
        self.simpli.fiber_bundles = self.fiber_bundles
        self.simpli.fiber_bundles_properties = self.fiber_bundles_properties

    def test_dict(self):
        d = self.simpli.get_dict()
        self.assertWarns(UserWarning, self.simpli.set_dict, d)

        # io
        with h5py.File('/tmp/fastpli.test.h5', 'w') as h5f:
            h5f['dict'] = str(d)
            d_read = h5f['dict'][...]
            d_new = dict(eval(str(d_read)))
            self.assertWarns(UserWarning, self.simpli.set_dict, d_new)

        self.addCleanup(os.remove, '/tmp/fastpli.test.h5')

    def test_dimension(self):
        self.simpli.fiber_bundles = [[[[1, 3, 0, 2], [1, 3, 7, 2]]]]
        self.simpli.fiber_bundles_properties = [[(1, 0, 0, 'p')]]
        self.simpli.dim = [3, 5, 7]
        self.simpli.dim_origin = [0, 0, 0]
        self.simpli.voxel_size = 1.0

        label_field, vec_field, tissue_properties = self.simpli.generate_tissue(
        )

        self.assertTrue(np.array_equal(label_field.shape, [3, 5, 7]))
        self.assertTrue(np.array_equal(vec_field.shape, [3, 5, 7, 3]))
        self.assertTrue(np.array_equal(tissue_properties.shape, [2, 2]))
        self.assertTrue(np.array_equal(label_field.shape, vec_field.shape[:3]))
        self.assertTrue(
            np.array_equal(label_field,
                           np.array(vec_field[:, :, :, 2], dtype=int)))

    def test_generator(self):
        self.simpli.dim = [10, 10, 10]
        self.simpli.voxel_size = 0.2
        label_field_0, vec_field, tissue_properties = self.simpli.generate_tissue(
            only_label=True)
        self.assertTrue(vec_field.size == 0)

        label_field_1, vec_field, tissue_properties = self.simpli.generate_tissue(
            only_label=False)
        self.assertTrue(np.array_equal(tissue_properties.shape, [4, 2]))
        self.assertTrue(np.array_equal(label_field_0, label_field_1))
        self.assertTrue(np.array_equal(label_field_1.shape,
                                       vec_field.shape[:3]))

    def test_cell_population(self):
        self.simpli.voxel_size = 1
        self.simpli.dim = [10, 10, 10]
        self.simpli.dim_origin = self.simpli.dim / 2
        self.simpli.cells_populations = [[[[0, 0, 0, 100]]]]
        self.simpli.cells_populations_properties = [[1, 10]]

        self.simpli.fiber_bundles = None
        self.simpli.fiber_bundles_properties = None

        label_field, vec_field, tissue_properties = self.simpli.generate_tissue(
        )

        self.assertTrue(np.all(label_field == 1))
        self.assertTrue(np.all(vec_field == 0))
        self.assertTrue(np.array_equal(tissue_properties.shape, [2, 2]))

    def test_simulator(self):
        self.simpli.voxel_size = 1
        self.simpli.dim = [10, 10, 10]
        self.simpli.dim_origin = self.simpli.dim / 2
        self.simpli.fiber_bundles = [[[[0, 0, 30, 100], [640, 640, 30, 100]]]]
        self.simpli.fiber_bundles_properties = [[(0.333, -0.004, 10, 'p'),
                                                 (0.666, -0.004, 5, 'b'),
                                                 (1.0, 0.004, 1, 'r')]]

        label_field, vec_field, tissue_properties = self.simpli.generate_tissue(
        )

        self.assertTrue(np.array_equal(tissue_properties.shape, [4, 2]))

        # PliSimulation ###
        self.simpli.filter_rotations = np.deg2rad([0, 30, 60, 90, 120, 150])
        self.simpli.light_intensity = 26000
        self.simpli.voxel_size = 1
        self.simpli.wavelength = 525
        self.simpli.untilt_sensor_view = False

        image = self.simpli.run_simulation(label_field, vec_field,
                                           tissue_properties, 0, 0)

        self.simpli.untilt_sensor_view = True
        image = self.simpli.run_simulation(label_field, vec_field,
                                           tissue_properties, 0, 0)

        with h5py.File('/tmp/fastpli.test.h5', 'w') as h5f:
            h5f['tissue'] = label_field.astype(np.uint16)
            h5f['vectorfield'] = vec_field
            h5f['data/0'] = image

        self.addCleanup(os.remove, '/tmp/fastpli.test.h5')

    def test_pipelline(self):
        self.simpli.voxel_size = 1
        self.simpli.dim = [10, 10, 10]
        self.simpli.dim_origin = self.simpli.dim / 2
        self.simpli.fiber_bundles = [[[[0, 0, 30, 100], [640, 640, 30, 100]]]]
        self.simpli.fiber_bundles_properties = [[(0.333, -0.004, 10, 'p'),
                                                 (0.666, -0.004, 5, 'b'),
                                                 (1.0, 0.004, 1, 'r')]]
        self.simpli.filter_rotations = np.deg2rad([0, 30, 60, 90, 120, 150])
        self.simpli.light_intensity = 26000
        self.simpli.voxel_size = 1
        self.simpli.wavelength = 525
        self.simpli.untilt_sensor_view = False
        self.simpli.sensor_gain = 3
        self.simpli.optical_sigma = 0.71
        self.simpli.resolution = 2
        self.simpli.tilts = np.deg2rad([(0, 0), (5.5, 0), (5.5, 90), (5.5, 180),
                                        (5.5, 270)])

        with h5py.File('/tmp/fastpli.test.h5', 'w') as h5f:
            with open(os.path.abspath(__file__), 'r') as script:
                self.simpli.run_pipeline(h5f=h5f,
                                         script=script.read(),
                                         save=["label_field", "vector_field"])

        with h5py.File('/tmp/fastpli.test.h5', 'r') as h5f:
            parameters = dict(eval(str(h5f['parameter/dict'][...])))
            self.assertTrue(self.simpli.get_dict() == parameters)
            self.assertWarns(UserWarning, self.simpli.set_dict, parameters)

        self.addCleanup(os.remove, '/tmp/fastpli.test.h5')

    def test_rofl(self):
        with self.assertRaisesRegex(ValueError, "tilts not set"):
            self.simpli.apply_rofl(0)

        with self.assertRaisesRegex(ValueError, "tilts not suitable for ROFL"):
            self.simpli.tilts = [0, 1]
            self.simpli.apply_rofl(0)

    def test_crop_tilt(self):
        self.simpli.voxel_size = 0.51
        self.simpli.resolution = 10
        self.simpli.set_voi([-100] * 3, [100] * 3)
        self.simpli.tilts = np.deg2rad([(15.5, 0)])
        dim_org = self.simpli.dim

        self.simpli.add_crop_tilt_halo()
        dim_crop = self.simpli.dim
        dim_crop[:2] -= 2 * self.simpli.crop_tilt_voxel()

        self.assertTrue(np.array_equal(dim_org, dim_crop))


if __name__ == '__main__':
    unittest.main()
