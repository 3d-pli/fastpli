import unittest
import numpy as np
import subprocess
import warnings
import h5py
import sys
import os
import tempfile

from fastpli.simulation import Simpli


class MainTest(unittest.TestCase):

    def setUp(self):
        self.fiber_bundles = [[[[0, 0, 0, 1], [1, 1, 1, 1], [2, 2, 2, 1]]]]
        self.fiber_bundles_properties = [[(0.333, -0.004, 10, 'p'),
                                          (0.666, 0, 5, 'b'),
                                          (1.0, 0.004, 1, 'r')]]

        self.simpli = Simpli()
        self.simpli.fiber_bundles = self.fiber_bundles
        self.simpli.fiber_bundles_properties = self.fiber_bundles_properties

    def test_dict(self):
        d = self.simpli.get_dict()
        self.assertWarns(UserWarning, self.simpli.set_dict, d)

        # # io
        # with h5py.File('/tmp/fastpli.test.h5', 'w') as h5f:
        #     h5f['dict'] = str(d)
        #     d_read = h5f['dict'][...]
        #     d_new = dict(eval(str(d_read)))
        #     self.assertWarns(UserWarning, self.simpli.set_dict, d_new)

        # self.addCleanup(os.remove, '/tmp/fastpli.test.h5')

    def test_dimension(self):
        self.simpli.fiber_bundles = [[[[1, 3, 0, 2], [1, 3, 7, 2]]]]
        self.simpli.fiber_bundles_properties = [[(1, 0, 0, 'p')]]
        self.simpli.dim = [3, 5, 7]
        self.simpli.dim_origin = [0, 0, 0]
        self.simpli.voxel_size = 1.0

        tissue, optical_axis, tissue_properties = self.simpli.generate_tissue()

        self.assertTrue(np.array_equal(tissue.shape, [3, 5, 7]))
        self.assertTrue(np.array_equal(optical_axis.shape, [3, 5, 7, 3]))
        self.assertTrue(np.array_equal(tissue_properties.shape, [2, 2]))
        self.assertTrue(np.array_equal(tissue.shape, optical_axis.shape[:3]))
        self.assertTrue(
            np.array_equal(tissue, np.array(optical_axis[:, :, :, 2],
                                            dtype=int)))

    def test_set_get_voi(self):
        with self.assertRaisesRegex(TypeError, "voxel_size is not set yet"):
            self.simpli.set_voi([0, 0, 0], [12, 12, 12])

        # TODO:
        pass

    def test_generator(self):
        self.simpli.dim = [10, 10, 10]
        self.simpli.voxel_size = 0.2
        tissue_0, optical_axis, tissue_properties = self.simpli.generate_tissue(
            only_tissue=True)
        self.assertTrue(optical_axis.size == 0)

        tissue_1, optical_axis, tissue_properties = self.simpli.generate_tissue(
            only_tissue=False)
        self.assertTrue(np.array_equal(tissue_properties.shape, [4, 2]))
        self.assertTrue(np.array_equal(tissue_0, tissue_1))
        self.assertTrue(np.array_equal(tissue_1.shape, optical_axis.shape[:3]))

    def test_tissue_radial(self):
        self.fiber_bundles = [[[[0, 0.25, -500, 1], [0, 0.25, 500, 1]]]]
        self.fiber_bundles_properties = [[(1.0, 0.004, 1, 'r')]]

        self.simpli = Simpli()
        self.simpli.fiber_bundles = self.fiber_bundles
        self.simpli.fiber_bundles_properties = self.fiber_bundles_properties

        self.simpli.voxel_size = 0.1
        self.simpli.set_voi([0, 0, 0], [1.5, 1.5, 1.5])

        tissue, optical_axis, tissue_properties = self.simpli.generate_tissue()

        self.assertTrue(tissue.dtype == np.int32)
        self.assertTrue(optical_axis.dtype == np.float32)

        for i in range(tissue.shape[-1]):
            self.assertTrue(np.array_equal(tissue[:, :, 0], tissue[:, :, i]))
            np.testing.assert_array_almost_equal(
                optical_axis[:, :, 0, :],
                optical_axis[:, :, i, :],
                decimal=np.finfo(optical_axis.dtype).precision)

        x = self.fiber_bundles[0][0][0, 0] / self.simpli.voxel_size
        y = self.fiber_bundles[0][0][0, 1] / self.simpli.voxel_size
        r = self.fiber_bundles[0][0][0, 3] / self.simpli.voxel_size

        for i in range(tissue.shape[0]):
            for j in range(tissue.shape[1]):
                t = 0 if (x - (i + 0.5))**2 + (y - (j + 0.5))**2 > r**2 else 1
                self.assertTrue(tissue[i, j, 0] == t)

                if t:
                    phi = np.arctan2((j + 0.5) - y, (i + 0.5) - x)
                    np.testing.assert_array_almost_equal(
                        optical_axis[i, j, 0, :],
                        [np.cos(phi), np.sin(phi), 0],
                        decimal=np.finfo(optical_axis.dtype).precision,
                        err_msg=f"x:{x}, y:{y}, i:{i},j:{j}")

    def test_tissue_parallel(self):
        self.fiber_bundles = [[[[0, 0.25, -500, 1], [0, 0.25, 500, 1]]]]
        self.fiber_bundles_properties = [[(1.0, -0.004, 1, 'p')]]

        self.simpli = Simpli()
        self.simpli.fiber_bundles = self.fiber_bundles
        self.simpli.fiber_bundles_properties = self.fiber_bundles_properties

        self.simpli.voxel_size = 0.1
        self.simpli.set_voi([0, 0, 0], [1.5, 1.5, 1.5])

        tissue, optical_axis, tissue_properties = self.simpli.generate_tissue()

        self.assertTrue(tissue.dtype == np.int32)
        self.assertTrue(optical_axis.dtype == np.float32)

        for i in range(tissue.shape[-1]):
            self.assertTrue(np.array_equal(tissue[:, :, 0], tissue[:, :, i]))
            np.testing.assert_array_almost_equal(
                optical_axis[:, :, 0, :],
                optical_axis[:, :, i, :],
                decimal=np.finfo(optical_axis.dtype).precision)

        x = self.fiber_bundles[0][0][0, 0] / self.simpli.voxel_size
        y = self.fiber_bundles[0][0][0, 1] / self.simpli.voxel_size
        r = self.fiber_bundles[0][0][0, 3] / self.simpli.voxel_size

        for i in range(tissue.shape[0]):
            for j in range(tissue.shape[1]):
                t = 0 if (x - (i + 0.5))**2 + (y - (j + 0.5))**2 > r**2 else 1
                self.assertTrue(tissue[i, j, 0] == t)

                if t:
                    np.testing.assert_array_almost_equal(
                        optical_axis[i, j, 0, :], [0, 0, 1],
                        decimal=np.finfo(optical_axis.dtype).precision,
                        err_msg=f"x:{x}, y:{y}, i:{i},j:{j}")

    def test_tissue_layered(self):
        self.fiber_bundles = [[[[0, 0.25, -500, 1], [0, 0.25, 500, 1]]]]
        self.fiber_bundles_properties = [[(0.3, -0.004, 1, 'p'),
                                          (0.6, 0, 1, 'b'),
                                          (1.0, 0.004, 1, 'r')]]

        self.simpli = Simpli()
        self.simpli.fiber_bundles = self.fiber_bundles
        self.simpli.fiber_bundles_properties = self.fiber_bundles_properties

        self.simpli.voxel_size = 0.1
        self.simpli.set_voi([0, 0, 0], [1.5, 1.5, 1.5])

        tissue, optical_axis, tissue_properties = self.simpli.generate_tissue()

        self.assertTrue(tissue.dtype == np.int32)
        self.assertTrue(optical_axis.dtype == np.float32)

        for i in range(tissue.shape[-1]):
            self.assertTrue(np.array_equal(tissue[:, :, 0], tissue[:, :, i]))
            np.testing.assert_array_almost_equal(
                optical_axis[:, :, 0, :],
                optical_axis[:, :, i, :],
                decimal=np.finfo(optical_axis.dtype).precision)

        x = self.fiber_bundles[0][0][0, 0] / self.simpli.voxel_size
        y = self.fiber_bundles[0][0][0, 1] / self.simpli.voxel_size
        r = self.fiber_bundles[0][0][0, 3] / self.simpli.voxel_size

        for i in range(tissue.shape[0]):
            for j in range(tissue.shape[1]):
                rr = (x - (i + 0.5))**2 + (y - (j + 0.5))**2
                t = 0 if rr > r**2 else 1

                if t:
                    if np.sqrt(rr) < 0.3 * r:
                        self.assertTrue(tissue[i, j, 0] == 1)
                        np.testing.assert_array_almost_equal(
                            optical_axis[i, j, 0, :], [0, 0, 1],
                            decimal=np.finfo(optical_axis.dtype).precision,
                            err_msg=f"x:{x}, y:{y}, i:{i},j:{j}")
                    elif np.sqrt(rr) < 0.6 * r:
                        self.assertTrue(tissue[i, j, 0] == 2)
                        self.assertTrue(
                            np.array_equal(optical_axis[i, j, 0, :], [0, 0, 0]))
                    else:
                        self.assertTrue(tissue[i, j, 0] == 3)
                        phi = np.arctan2((j + 0.5) - y, (i + 0.5) - x)
                        np.testing.assert_array_almost_equal(
                            optical_axis[i, j, 0, :],
                            [np.cos(phi), np.sin(phi), 0],
                            decimal=np.finfo(optical_axis.dtype).precision,
                            err_msg=f"x:{x}, y:{y}, i:{i},j:{j}")
                else:
                    self.assertTrue(tissue[i, j, 0] == t)

        self.simpli.voxel_size = 1
        self.simpli.set_voi([-10] * 3, [10] * 3)
        self.simpli.fiber_bundles = [[[[-100, 0, 0, 10], [100, 0, 0, 10]]]]
        self.simpli.fiber_bundles_properties = [[(0.333, -0.004, 10, 'p'),
                                                 (0.666, 0, 5, 'b'),
                                                 (1.0, 0.004, 1, 'r')]]

        tissue, optical_axis, tissue_properties = self.simpli.generate_tissue()

        self.assertTrue(np.any(tissue != 0))
        for x in range(tissue.shape[0]):
            for y in range(tissue.shape[0]):
                for z in range(tissue.shape[0]):
                    if tissue[x, y, z] == 2:
                        self.assertTrue(np.all(optical_axis[x, y, z, :] == 0))
                    elif tissue[x, y, z]:
                        self.assertTrue(np.any(optical_axis[x, y, z, :] != 0))
                    else:
                        self.assertTrue(np.all(optical_axis[x, y, z, :] == 0))

        self.assertTrue(np.array_equal(tissue_properties.shape, [4, 2]))

    def test_cell_population(self):
        self.simpli.voxel_size = 1
        self.simpli.dim = [10, 10, 10]
        self.simpli.dim_origin = self.simpli.dim / 2
        self.simpli.cells_populations = [[[[0, 0, 0, 100]]]]
        self.simpli.cells_populations_properties = [[1, 10]]

        self.simpli.fiber_bundles = None
        self.simpli.fiber_bundles_properties = None

        tissue, optical_axis, tissue_properties = self.simpli.generate_tissue()

        self.assertTrue(np.all(tissue == 1))
        self.assertTrue(np.all(optical_axis == 0))
        self.assertTrue(np.array_equal(tissue_properties.shape, [2, 2]))

    def test_simulator(self):
        self.simpli.voxel_size = 1
        self.simpli.set_voi([-10] * 3, [10] * 3)
        self.simpli.fiber_bundles = [[[[-100, 0, 0, 10], [100, 0, 0, 10]]]]
        self.simpli.fiber_bundles_properties = [[(0.333, -0.004, 10, 'p'),
                                                 (0.666, 0, 5, 'b'),
                                                 (1.0, 0.004, 1, 'r')]]

        tissue, optical_axis, tissue_properties = self.simpli.generate_tissue()

        # PliSimulation ###
        self.simpli.filter_rotations = np.deg2rad([0, 30, 60, 90, 120, 150])
        self.simpli.light_intensity = 26000
        self.simpli.voxel_size = 1
        self.simpli.wavelength = 525

        # normal view
        self.simpli.untilt_sensor_view = False
        image = self.simpli.run_simulation(tissue, optical_axis,
                                           tissue_properties, 0, 0)
        self.assertFalse(np.any(np.isnan(image)))
        self.assertTrue(np.any(image != self.simpli.light_intensity / 4))

        image = self.simpli.run_simulation(tissue, optical_axis,
                                           tissue_properties, np.deg2rad(8),
                                           np.deg2rad(42))
        self.assertFalse(np.any(np.isnan(image)))
        self.assertTrue(np.any(image != self.simpli.light_intensity / 4))

        # untilt view
        self.simpli.untilt_sensor_view = True
        image = self.simpli.run_simulation(tissue, optical_axis,
                                           tissue_properties, 0, 0)
        self.assertFalse(np.any(np.isnan(image)))
        self.assertTrue(np.any(image != self.simpli.light_intensity / 4))

        image = self.simpli.run_simulation(tissue, optical_axis,
                                           tissue_properties, np.deg2rad(8),
                                           np.deg2rad(42))
        self.assertFalse(np.any(np.isnan(image)))
        self.assertTrue(np.any(image != self.simpli.light_intensity / 4))

    def test_pipeline(self):
        self.simpli.voxel_size = 1
        self.simpli.dim = [10, 10, 10]
        self.simpli.dim_origin = self.simpli.dim / 2
        self.simpli.fiber_bundles = [[[[0, 0, 30, 100], [640, 640, 30, 100]]]]
        self.simpli.fiber_bundles_properties = [[(0.333, -0.004, 10, 'p'),
                                                 (0.666, 0, 5, 'b'),
                                                 (1.0, 0.004, 1, 'r')]]
        self.simpli.filter_rotations = np.deg2rad([0, 30, 60, 90, 120, 150])
        self.simpli.light_intensity = 26000
        self.simpli.voxel_size = 1
        self.simpli.wavelength = 525
        self.simpli.untilt_sensor_view = False
        self.simpli.optical_sigma = 0.71
        self.simpli.noise_model = lambda x: np.random.negative_binomial(
            x / (3 - 1), 1 / 3)
        self.simpli.pixel_size = 2
        self.simpli.tilts = np.deg2rad([(0, 0), (5.5, 0), (5.5, 90), (5.5, 180),
                                        (5.5, 270)])

        self.simpli.run_pipeline(save=["tissue", "optical_axis"])

        with h5py.File(os.path.join(tempfile.gettempdir(), 'fastpli.test.h5'),
                       'w') as h5f:
            with open(os.path.abspath(__file__), 'r') as script:
                self.simpli.run_pipeline(h5f=h5f,
                                         script=script.read(),
                                         save=["tissue", "optical_axis"])

        # with h5py.File('/tmp/fastpli.test.h5', 'r') as h5f:
        #     parameters = dict(eval(str(h5f.attrs['fastpli/simpli'])))
        #     self.assertTrue(self.simpli.get_dict() == parameters)
        #     self.assertWarns(UserWarning, self.simpli.set_dict, parameters)

        self.addCleanup(os.remove,
                        os.path.join(tempfile.gettempdir(), 'fastpli.test.h5'))

    def test_rofl(self):
        with self.assertRaisesRegex(ValueError, "tilts not set"):
            self.simpli.apply_rofl(0)

        with self.assertRaisesRegex(ValueError, "tilts not suitable for ROFL"):
            self.simpli.tilts = [0, 1]
            self.simpli.apply_rofl(0)

    def test_crop_tilt(self):
        self.simpli.voxel_size = 0.51
        self.simpli.pixel_size = 10
        self.simpli.set_voi([-100] * 3, [100] * 3)
        self.simpli.tilts = np.deg2rad([(15.5, 0)])
        dim_org = self.simpli.dim.copy()

        self.simpli.add_crop_tilt_halo()
        dim_crop = self.simpli.dim.copy()
        dim_crop[:2] -= 2 * self.simpli.crop_tilt_voxel()
        self.assertTrue(np.array_equal(dim_org, dim_crop))

        image_org = np.empty(dim_org)
        image_halo = np.empty(self.simpli.dim)
        image_crop = self.simpli.rm_crop_tilt_halo(image_halo)
        self.assertTrue(np.array_equal(image_org.shape, image_crop.shape))

        self.simpli.tilts = np.deg2rad([(0, 100)])
        dim_org = self.simpli.dim.copy()
        self.simpli.add_crop_tilt_halo()
        dim_crop = self.simpli.dim.copy()
        self.assertTrue(np.array_equal(dim_org, dim_crop))
        image_org = np.empty(dim_org)
        image_halo = np.empty(self.simpli.dim)
        image_crop = self.simpli.rm_crop_tilt_halo(image_halo)
        self.assertTrue(np.array_equal(image_org.shape, image_crop.shape))

    def test_length_vec(self):

        self.simpli.voxel_size = 1
        self.simpli.set_voi([-10] * 3, [10] * 3)
        self.simpli.fiber_bundles = [[[[-100, 0, 0, 10], [100, 0, 0, 10]]]]
        self.simpli.fiber_bundles_properties = [[(0.333, -0.004, 10, 'p'),
                                                 (0.666, 0, 5, 'b'),
                                                 (1.0, 0.004, 1, 'r')]]

        tissue, optical_axis, tissue_properties = self.simpli.generate_tissue()

        # PliSimulation ###
        self.simpli.filter_rotations = np.deg2rad([0, 30, 60, 90, 120, 150])
        self.simpli.light_intensity = 26000
        self.simpli.voxel_size = 1
        self.simpli.wavelength = 525

        # normal view
        self.simpli.untilt_sensor_view = False

        with warnings.catch_warnings(record=True) as w:
            self.simpli.run_simulation(tissue, optical_axis, tissue_properties,
                                       0, 0)
            self.assertTrue(len(w) == 0)

        optical_axis[tissue > 0, :] *= 1.0001
        self.assertWarns(UserWarning, self.simpli.run_simulation, tissue,
                         optical_axis, tissue_properties, 0, 0)

        optical_axis[:] = 0
        self.assertWarns(UserWarning, self.simpli.run_simulation, tissue,
                         optical_axis, tissue_properties, 0, 0)

    def test_prev_result(self):
        FILE_NAME = os.path.abspath(__file__)
        FILE_PATH = os.path.dirname(FILE_NAME)
        subprocess.run(
            [sys.executable,
             os.path.join(FILE_PATH, "simpli_rep.py")],
            stdout=subprocess.DEVNULL,
            check=True)
        result = subprocess.run([
            "h5diff",
            "--relative=0.00001",  # some have 32bit
            os.path.join(FILE_PATH, "simpli_rep.h5"),
            os.path.join(tempfile.gettempdir(), 'simpli_rep.h5')
        ]).returncode == 0
        if not result:
            subprocess.run([
                "h5diff", "--relative=0.00001", "-r",
                os.path.join(FILE_PATH, "simpli_rep.h5"),
                os.path.join(tempfile.gettempdir(), 'simpli_rep.h5')
            ])
        self.assertTrue(result)

        self.addCleanup(
            os.remove,
            os.path.join(os.path.join(tempfile.gettempdir(), "simpli_rep.h5")))


if __name__ == '__main__':
    unittest.main()
