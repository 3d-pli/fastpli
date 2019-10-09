# from . import __generation as generation
from .__generation import _Generator
from .__simulation import _Simulator

from ..version import __version__

from fastpli.analysis import rofl, epa, affine_transformation
from fastpli.simulation import optic
from fastpli.tools import rotation

import numpy as np
import warnings
import copy

import h5py
from PIL import Image

# TODO: write json -> parameter function


class Simpli:
    __isfrozen = False

    def __setattr__(self, key, value):
        if self.__isfrozen and not hasattr(self, key):
            raise TypeError("%r is a frozen class" % self)
        object.__setattr__(self, key, value)

    def _freeze(self):
        self.__isfrozen = True

    def __init__(self, mpi_comm=None):

        self._gen = _Generator()
        self._sim = _Simulator()
        if mpi_comm:
            from mpi4py import MPI
            self._gen.set_mpi_comm(MPI._addressof(mpi_comm))
            self._sim.set_mpi_comm(MPI._addressof(mpi_comm))

        self._fiber_bundles = None
        self._fiber_bundles_properties = None
        self._cells_populations = None
        self._cells_populations_properties = None
        self._dim = None
        self._dim_origin = np.array([0, 0, 0], dtype=np.float64)
        self._voxel_size = None
        self._resolution = None

        self._filter_rotations = None
        self._light_intensity = None
        self._step_size = 1.0
        self._interpolate = True
        self._untilt_sensor_view = True
        self._flip_z_beam = False
        self._wavelength = None
        self._tissue_refrection = 1

        self._omp_num_threads = 1
        self._debug = False

        self._freeze()

    def as_dict(self):
        return {
            'version': __version__,
            # 'fiber_bundles': self._fiber_bundles,
            'fiber_bundles_properties': self._fiber_bundles_properties,
            # 'cells_populations ': self._cells_populations,
            'cells_populations_properties ': self._cells_populations_properties,
            'dim': self._dim,
            'dim_origin': self._dim_origin,
            'voxel_size': self._voxel_size,
            'resolution': self._resolution,
            'filter_rotations': self._filter_rotations,
            'light_intensity': self._light_intensity,
            'step_size': self._step_size,
            'interpolate': self._interpolate,
            'untilt_sensor_view': self._untilt_sensor_view,
            'flip_z_beam': self._flip_z_beam,
            'wavelength': self._wavelength,
            'tissue_refrection': self._tissue_refrection,
            'omp_num_threads': self._omp_num_threads,
            'debug': self._debug
        }

    def _print_debug(self, msg):
        if self._debug:
            print("debug: " + msg)

    @property
    def debug(self):
        return self._debug

    @debug.setter
    def debug(self, debug):
        self._debug = bool(debug)

    @property
    def dim(self):
        return self._dim

    @dim.setter
    def dim(self, dim):
        self._dim = np.array(dim, dtype=np.int64)

        if not np.array_equal(self._dim, dim):
            raise ValueError("dim would changed because it was not an integer")

    @property
    def dim_origin(self):
        return self._dim_origin

    @dim_origin.setter
    def dim_origin(self, dim_origin):
        self._dim_origin = np.array(dim_origin, dtype=np.float64)

    @property
    def voxel_size(self):
        return (self._voxel_size)

    @voxel_size.setter
    def voxel_size(self, voxel_size):
        if not isinstance(voxel_size, (int, float)):
            raise TypeError("voxel_size : (int, float)")

        if voxel_size <= 0:
            raise ValueError("voxel_size <= 0")

        self._voxel_size = float(voxel_size)

    @property
    def resolution(self):
        return (self._resolution)

    @resolution.setter
    def resolution(self, resolution):
        if not isinstance(resolution, (int, float)):
            raise TypeError("resolution : (int, float)")
        if resolution <= 0:
            raise ValueError("resolution <= 0")

        self._resolution = float(resolution)

    @property
    def voi(self):
        raise NotImplementedError("voi has to be get/set via get/set_voi()")

    @voi.setter
    def voi(self, voi):
        raise NotImplementedError("voi has to be get/set via get/set_voi()")

    def get_voi(self):
        if self._voxel_size is None:
            return None
            self._print_debug("voxel_size is not set, voi can't be calculated")

        if self._dim is None:
            return None
            self._print_debug("dim is not set, voi can't be calculated")

        if self._dim_origin is None:
            return None
            self._print_debug("dim_origin is not set, voi can't be calculated")

        voi = np.zeros((6,))
        voi[::2] = self._dim_origin
        voi[1::2] = voi[::2] + self._dim * self._voxel_size
        return voi

    def set_voi(self, voi):
        voi = np.array(voi, dtype=np.float64)
        if voi.size != 6 or voi.shape[0] != 6:
            raise TypeError("voi: wrong shape, has to be (6,)")

        if voi[0] > voi[1] or voi[2] > voi[3] or voi[4] > voi[5]:
            raise ValueError(
                "voi not corrected sorted: (x_min, x_max, y_min, y_max, z_min, z_max)"
            )

        if self._voxel_size is None:
            raise TypeError("voxel_size is not set yet")
        tmp = np.array(voi / self._voxel_size)
        self._dim = np.array(
            (int(tmp[1] - tmp[0]), int(tmp[3] - tmp[2]), int(tmp[5] - tmp[4])),
            dtype=np.int64)
        self._dim_origin = voi[::2]

        self._print_debug("dim and dim_origin recalculated")

    @property
    def filter_rotations(self):
        return (self._filter_rotations)

    @filter_rotations.setter
    def filter_rotations(self, filter_rotations):

        filter_rotations = np.array(filter_rotations, dtype=np.float64)

        if filter_rotations.size == 0 or filter_rotations.ndim != 1:
            raise TypeError("filter_rotations : nx1")

        self._filter_rotations = filter_rotations

    @property
    def light_intensity(self):
        return self._light_intensity

    @light_intensity.setter
    def light_intensity(self, light_intensity):
        self._light_intensity = float(light_intensity)

    @property
    def wavelength(self):
        return self._wavelength

    @wavelength.setter
    def wavelength(self, wavelength):
        self._wavelength = float(wavelength)

    @property
    def tissue_refrection(self):
        return self._tissue_refrection

    @tissue_refrection.setter
    def tissue_refrection(self, tissue_refrection):
        self._tissue_refrection = float(tissue_refrection)

    @property
    def step_size(self):
        return self._step_size

    @step_size.setter
    def step_size(self, step_size):
        self._step_size = float(step_size)

    @property
    def interpolate(self):
        return self._interpolate

    @interpolate.setter
    def interpolate(self, interpolate):
        self._interpolate = bool(interpolate)

    @property
    def untilt_sensor_view(self):
        return self._untilt_sensor_view

    @untilt_sensor_view.setter
    def untilt_sensor_view(self, untilt_sensor_view):
        self._untilt_sensor_view = bool(untilt_sensor_view)

    @property
    def flip_z_beam(self):
        return self._flip_z_beam

    @flip_z_beam.setter
    def flip_z_beam(self, flip_z_beam):
        self._flip_z_beam = bool(flip_z_beam)

    @property
    def fiber_bundles(self):
        return self._fiber_bundles

    @fiber_bundles.setter
    def fiber_bundles(self, fbs):

        if fbs is None:
            self._fiber_bundles = None
            return

        if not isinstance(fbs, (list, tuple)):
            raise TypeError("fbs is not a list")

        for fb_i, fb in enumerate(fbs):
            if not isinstance(fb, (list, tuple)):
                raise TypeError("fb is not a list")

            for f_i, f in enumerate(fb):
                fbs[fb_i][f_i] = np.array(f, dtype=np.float64, copy=False)

                if len(fbs[fb_i][f_i].shape) != 2:
                    raise TypeError("fiber size need to be nx4")

                if fbs[fb_i][f_i].shape[1] != 4:
                    raise TypeError("fiber size need to be nx4")

        self._fiber_bundles = fbs

    @property
    def fiber_bundles_properties(self):
        return self._fiber_bundles_properties

    @fiber_bundles_properties.setter
    def fiber_bundles_properties(self, bundle_layer_properties):

        if bundle_layer_properties is None:
            self._fiber_bundles_properties = None
            return

        if not isinstance(bundle_layer_properties, (list, tuple)):
            raise TypeError("properties must be a list(list(tuples))")

        if not self._fiber_bundles:
            raise ValueError("fiber_bundles have not been set yet")

        if len(self._fiber_bundles) != len(bundle_layer_properties):
            raise TypeError(
                "properties must have the same size as fiber_bundles")

        self._fiber_bundles_properties = []
        for prop in bundle_layer_properties:
            if not isinstance(prop, (list, tuple)):
                raise TypeError("properties must be a list(list(tuples))")

            self._fiber_bundles_properties.append([])
            for ly in prop:
                if len(ly) != 4:
                    raise TypeError(
                        "properties must have len 4 (float, float, float, char)"
                    )

                if ly[1] < 0 and ly[-1] is 'r':
                    warnings.warn('birefringence negative and radial')
                if ly[1] > 0 and ly[-1] is 'p':
                    warnings.warn('birefringence positive and parallel')

        self._fiber_bundles_properties = bundle_layer_properties

    @property
    def cells_populations(self):
        return self._cells_populations

    @cells_populations.setter
    def cells_populations(self, cps):
        if cps is None:
            self._cells_populations = None
            return

        if not isinstance(cps, (list, tuple)):
            raise TypeError("cells_populations is not a list")

        for cp_i, cp in enumerate(cps):
            if not isinstance(cp, (list, tuple)):
                raise TypeError("cells_population is not a list")

            for c_i, c in enumerate(cp):
                cps[cp_i][c_i] = np.array(c, dtype=np.float64, copy=False)

                if len(cps[cp_i][c_i].shape) != 2:
                    raise TypeError("cell size need to be nx4")

                if cps[cp_i][c_i].shape[1] != 4:
                    raise TypeError("cell size need to be nx4")

        self._cells_populations = cps

    @property
    def cells_populations_properties(self):
        return self._cells_populations_properties

    @cells_populations_properties.setter
    def cells_populations_properties(self, cells_populations_properties):

        if cells_populations_properties is None:
            self._cells_populations_properties = None
            return

        if not isinstance(cells_populations_properties, (list, tuple)):
            raise TypeError("properties must be a list")

        self._cells_populations_properties = []
        for prop in cells_populations_properties:

            if not isinstance(prop, (list, tuple)):
                raise TypeError("properties must be a list of 2 arguments")

            if len(prop) != 2:
                raise TypeError("properties must be a list of 2 arguments")

        self._cells_populations_properties = cells_populations_properties

    def _check_property_length(self):
        if self._fiber_bundles:
            if len(self._fiber_bundles) != len(self._fiber_bundles_properties):
                raise TypeError(
                    "properties must have the same size as fiber_bundles")

        if self._cells_populations:
            if len(self._cells_populations) != len(
                    self._cells_populations_properties):
                raise TypeError(
                    "properties must have the same size as cell_populations")

    def generate_tissue(self, only_label=False, progress_bar=False):

        if self._dim is None:
            raise ValueError('dim not set')

        if self._dim_origin is None:
            raise ValueError('dim_origin not set')

        if self._voxel_size is None:
            raise ValueError('voxel_size not set')

        self._check_property_length()

        self._gen.set_volume(self._dim, self._dim_origin, self._voxel_size)
        if self._fiber_bundles:
            self._gen.set_fiber_bundles(self._fiber_bundles,
                                        self._fiber_bundles_properties)
        if self._cells_populations:
            self._gen.set_cell_populations(self._cells_populations,
                                           self._cells_populations_properties)
        label_field, vec_field, tissue_properties = self._gen.run_generation(
            only_label, progress_bar)

        return label_field, vec_field, tissue_properties

    def _init_pli_setup(self):
        if self._step_size <= 0:
            raise ValueError('step_size <= 0')

        if self._light_intensity is None:
            raise ValueError('light_intensity not set')

        if self._voxel_size is None:
            raise ValueError('voxel_size not set')

        if self._wavelength is None:
            raise ValueError('wavelength not set')

        if self._filter_rotations is None:
            raise ValueError('filter_rotations not set')

        self._sim.set_pli_setup(self._step_size, self._light_intensity,
                                self._voxel_size, self._wavelength,
                                self._tissue_refrection, self._interpolate,
                                self._untilt_sensor_view, self._flip_z_beam,
                                self._filter_rotations)

    def run_simulation(self, label_field, vec_field, tissue_properties, theta,
                       phi):

        label_field = np.array(label_field, dtype=np.int32, copy=False)
        vec_field = np.array(vec_field, dtype=np.float32, copy=False)

        self._init_pli_setup()

        tissue_properties = np.array(tissue_properties)
        if len(tissue_properties.shape) != 2:
            raise ValueError("tissue_properties wrong shape")

        if tissue_properties.shape[1] != 2:
            raise ValueError("tissue_properties wrong shape")

        if tissue_properties.shape[0] <= np.max(label_field.flatten()):
            raise ValueError(
                "tissue_properties.shape[0] < np.max(label_field.flatten())")

        image = self._sim.run_simulation(self._dim, label_field, vec_field,
                                         tissue_properties, theta, phi)

        return image.astype(
            np.float32)  # optic resize will force float32 because of PIL

    @property
    def omp_num_threads(self):
        return self._omp_num_threads

    @omp_num_threads.setter
    def omp_num_threads(self, num_threads):

        if not isinstance(num_threads, int):
            raise TypeError("num_threads has to be a positiv integer")

        if num_threads <= 0:
            raise TypeError("num_threads has to be a positiv integer")

        num_threads_gen = self._gen.set_omp_num_threads(num_threads)
        num_threads_sim = self._sim.set_omp_num_threads(num_threads)

        if num_threads_gen != num_threads_sim:
            raise ValueError("num_threads_gen != num_threads_sim")

        if num_threads_gen != num_threads:
            warnings.warn("reduced num_threads: " + str(num_threads_gen))

        self._omp_num_threads = num_threads_gen

    def memory_usage(self, unit='MB', item='all'):
        if not isinstance(item, str):
            raise TypeError('item has to be str')

        if not isinstance(unit, str):
            raise TypeError('unit has to be str')

        if unit == 'B':
            div = 1024**0
        elif unit == 'kB':
            div = 1024**1
        elif unit == 'MB':
            div = 1024**2
        elif unit == 'GB':
            div = 1024**3
        else:
            raise ValueError('allowed is only B, kB, MB, GB')

        if self._dim is None:
            raise TypeError('dimension not set yet')

        if item == 'label_field':
            return np.prod(self._dim) * (32 + 32) / 8 / div
        elif item == 'all':
            return np.prod(self._dim) * (32 + 32 + 3 * 32) / 8 / div
        else:
            raise ValueError('allowed is only label_field or all')

    def save_mpi_array_as_h5(self, h5f, data, data_name, lock_dim=None):

        dim_local = self._gen.dim_local()
        dim_offset = self._gen.dim_offset()

        if not isinstance(data, np.ndarray):
            raise TypeError(
                "only numpy arrays are compatible with save_mpi_array_as_h5")

        dset_dim = np.copy(self._dim)
        if len(data.shape) < len(dset_dim):
            dset_dim = dset_dim[:len(data.shape)]
        if len(data.shape) > len(dset_dim):
            dset_dim = np.append(dset_dim, data.shape[3:])

        if lock_dim:
            if isinstance(lock_dim, int):
                lock_dim = [lock_dim]

            lock_dim = list(lock_dim)
            for l in lock_dim:
                dset_dim[l] = data.shape[l]

        dset = h5f.create_dataset(data_name, dset_dim, dtype=data.dtype)

        if len(data.shape) == 2:
            if data.size * data.itemsize > 2 * (2**10)**3:  # 2 GB
                for i in range(data.shape[0]):
                    dset[i + dim_offset[0], dim_offset[1]:dim_offset[1] +
                         dim_local[1]] = data[i, :]
            else:
                dset[dim_offset[0]:dim_offset[0] +
                     dim_local[0], dim_offset[1]:dim_offset[1] +
                     dim_local[1]] = data

        elif len(data.shape) == 3:
            if data.size * data.itemsize > 2 * (2**10)**3:  # 2 GB
                for i in range(data.shape[0]):
                    dset[i + dim_offset[0], dim_offset[1]:dim_offset[1] +
                         dim_local[1], dim_offset[2]:dim_offset[2] +
                         dim_local[2]] = data[i, :]
            else:
                dset[dim_offset[0]:dim_offset[0] +
                     dim_local[0], dim_offset[1]:dim_offset[1] +
                     dim_local[1], dim_offset[2]:dim_offset[2] +
                     dim_local[2]] = data

        elif len(data.shape) > 3:
            if data.size * data.itemsize > 2 * (2**10)**3:  # 2 GB
                for i in range(data.shape[0]):
                    dset[i + dim_offset[0], dim_offset[1]:dim_offset[1] +
                         dim_local[1], dim_offset[2]:dim_offset[2] +
                         dim_local[2], :] = data[i, :]
            else:
                dset[dim_offset[0]:dim_offset[0] +
                     dim_local[0], dim_offset[1]:dim_offset[1] +
                     dim_local[1], dim_offset[2]:dim_offset[2] +
                     dim_local[2], :] = data

        else:
            raise TypeError("no compatible save_mpi_array_as_h5: " + data_name)

    def apply_optic(
            self,
            image_stack,
            delta_sigma=0.71,  # only for LAP!
            gain=3.0,  # only for LAP!
            resample_mode=Image.BILINEAR):

        if self._resolution is None:
            raise TypeError("resolution is not set")

        return optic.apply(image_stack, self._voxel_size, self._resolution,
                           delta_sigma, gain, resample_mode)

    def apply_resize(self, image, resample_mode=Image.BILINEAR):

        if self._resolution is None:
            raise TypeError("resolution is not set")

        if self._voxel_size is None:
            raise TypeError("voxel_size is not set")

        scale = self._resolution / self._voxel_size
        size = np.array(np.array(image.shape) // scale, dtype=int)
        return optic.resize(image, size, resample_mode)

    def apply_untilt(self, images, theta, phi, mode='nearest'):

        if theta == 0:
            return images

        # calculate transformation matrix
        p = self._dim
        p_rot = 0.5 * np.array([p[0], p[1], 0])
        p_out = np.array([[p[0], p[1], 0], [p[0], 0, 0], [0, p[1], 0]])
        rot = rotation.theta_phi(-theta, phi)

        p_in = np.array([np.dot(rot, p - p_rot) + p_rot for p in p_out])

        # TODO: refraction has to be implemented

        M = affine_transformation.calc_matrix(p_in[:, :2], p_out[:, :2])

        images_untilt = affine_transformation.image(images, M, mode)
        if images.ndim == 3:
            images_untilt = np.atleast_3d(images_untilt)

        return images_untilt

    def apply_epa(self, image_stack, mask=None):

        transmittance, direction, retardation = epa.epa(image_stack)
        if mask is not None:
            transmittance[np.invert(mask)] = float('nan')
            direction[np.invert(mask)] = float('nan')
            retardation[np.invert(mask)] = float('nan')

        return transmittance, direction, retardation

    def apply_rofl(
            self,
            tilting_stack,
            tilt_angle=np.deg2rad(5.5),  # only LAP!
            gain=3.0,  # only LAP!
            dir_offset=0,
            mask=None,
            num_threads=2,
            grad_mode=False):

        rofl_direction, rofl_incl, rofl_t_rel, dirdevmap, incldevmap, treldevmap, funcmap, itermap = rofl.rofl_stack(
            tilting_stack, tilt_angle, gain, dir_offset, mask, num_threads,
            grad_mode)

        return rofl_direction, rofl_incl, rofl_t_rel, (dirdevmap, incldevmap,
                                                       treldevmap, funcmap,
                                                       itermap)
