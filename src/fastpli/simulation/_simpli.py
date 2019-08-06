# from . import __generation as generation
from .__generation import _Generator, _LayerProperty
from .__simulation import _Simulator

from fastpli.analysis import rofl, epa
from fastpli.objects import Fiber, Cell
from fastpli.simulation import optic
from fastpli.tools import rotation

import numpy as np
import h5py
from PIL import Image

# TODO: write json -> parameter function


class Simpli:

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
        self._dim_origin = np.array([0, 0, 0], dtype=float)
        self._voxel_size = None
        self._resolution = None
        self._flip_direction = np.array([0, 0, 0], dtype=bool)

        self._filter_rotations = None
        self._light_intensity = None
        self._wavelength = None
        self._untilt_sensor = None

        self._debug = False

    def print_debug(self, msg):
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
        self._dim = np.array(dim, dtype=int)

        if not np.array_equal(self._dim, dim):
            raise ValueError("dim would changed because it was not an integer")

    @property
    def dim_origin(self):
        return self._dim_origin

    @dim_origin.setter
    def dim_origin(self, dim_origin):
        self._dim_origin = np.array(dim_origin, dtype=float)

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

    @voxel_size.setter
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
            self.print_debug("voxel_size is not set, voi can't be calculated")

        if self._dim is None:
            return None
            self.print_debug("dim is not set, voi can't be calculated")

        if self._dim_origin is None:
            return None
            self.print_debug("dim_origin is not set, voi can't be calculated")

        voi = np.zeros((6,))
        voi[::2] = self._dim_origin
        voi[1::2] = voi[::2] + self._dim * self._voxel_size
        return voi

    def set_voi(self, voi):
        voi = np.array(voi, dtype=float)
        if voi.size != 6 or voi.shape[0] != 6:
            raise TypeError("voi: wrong shape, has to be (6,)")

        if voi[0] > voi[1] or voi[2] > voi[3] or voi[4] > voi[5]:
            raise ValueError(
                "voi not corrected sorted: (x_min, x_max, y_min, y_max, z_min, z_max)"
            )

        self._voi = voi

        if self._voxel_size is None:
            self.print_debug(
                "voxel_size is not set yet, dim and dim_origin will be calculated when setted"
            )
            return

        tmp = np.array(self._voi / self._voxel_size)
        self._dim = np.array(
            (int(tmp[1] - tmp[0]), int(tmp[3] - tmp[2]), int(tmp[5] - tmp[4])),
            dtype=int)
        self._dim_origin = self._voi[::2]

        self.print_debug("dim and dim_origin recalculated")

    @property
    def flip_direction(self):
        return (self._flip_direction)

    @flip_direction.setter
    def flip_direction(self, flip_direction):

        flip_direction = np.array(flip_direction, dtype=bool)

        if flip_direction.size != 3 or flip_direction.shape != 1:
            raise TypeError("flip_direction : 3d bool")

        self._flip_direction = flip_direction

    @property
    def filter_rotations(self):
        return (self._filter_rotations)

    @filter_rotations.setter
    def filter_rotations(self, filter_rotations):

        filter_rotations = np.array(filter_rotations, dtype=float)

        if filter_rotations.size == 0 or filter_rotations.ndim != 1:
            print(filter_rotations)
            print(filter_rotations.shape)
            raise TypeError("filter_rotations : nx1")

        self._filter_rotations = filter_rotations

    @property
    def light_intensity(self):
        return self._light_intensity

    @light_intensity.setter
    def light_intensity(self, light_intensity):
        self._light_intensity = light_intensity

    @property
    def wavelength(self):
        return self._wavelength

    @wavelength.setter
    def wavelength(self, wavelength):
        self._wavelength = wavelength

    @property
    def untilt_sensor(self):
        return self._untilt_sensor

    @untilt_sensor.setter
    def untilt_sensor(self, untilt_sensor):
        self._untilt_sensor = untilt_sensor

    @property
    def fiber_bundles(self):
        return self._fiber_bundles

    @fiber_bundles.setter
    def fiber_bundles(self, fbs):
        if not isinstance(fbs, (list, tuple)):
            raise TypeError("fbs is not a list")

        for fb_i, fb in enumerate(fbs):
            if not isinstance(fb, (list, tuple)):
                raise TypeError("fb is not a list")

            for f_i, f in enumerate(fb):
                fbs[fb_i][f_i] = np.array(f, dtype=np.float32, copy=False)

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
                self._fiber_bundles_properties[-1].append(
                    _LayerProperty(ly[0], ly[1], ly[2], ly[3]))

    @property
    def cells_populations(self):
        return self._cells_populations

    @cells_populations.setter
    def cells_populations(self, cps):
        if not isinstance(cps, (list, tuple)):
            raise TypeError("cells_populations is not a list")

        for cp_i, cp in enumerate(cps):
            if not isinstance(cp, (list, tuple)):
                raise TypeError("cells_population is not a list")

            for c_i, c in enumerate(cp):
                if isinstance(c, Cell):
                    continue
                elif isinstance(c, (list, tuple)):
                    c = np.array(c)

                if isinstance(c, np.ndarray):
                    if len(c.shape) is not 2 or c.shape[1] is not 4:
                        raise TypeError("cell elements has to be of dim nx4")
                    cps[cp_i][c_i] = Cell(c[:, 0:-1], c[:, -1])
                else:
                    raise TypeError(
                        "cell hast to be a objects.Cell, 4d-list or 4d-array")

        self._cells_populations = cps

    @property
    def cells_populations_properties(self):
        return self._cells_populations_properties

    @cells_populations_properties.setter
    def cells_populations_properties(self, cells_populations_properties):

        if not isinstance(cells_populations_properties, (list, tuple)):
            raise TypeError("properties must be a list")

        self._cells_populations_properties = []
        for prop in cells_populations_properties:

            if not isinstance(prop, (list, tuple)):
                raise TypeError("properties must be a list of 2 arguments")

            if len(prop) != 2:
                raise TypeError("properties must be a list of 2 arguments")

            self._cells_populations_properties.append(
                generation.__CellProperty(prop[0], prop[1]))

    def ReadFiberFile(self, filename):
        self._fiber_bundles = []
        with h5py.File(filename, 'r') as h5f:

            fbs = h5f['fiber_bundles']
            for fb in fbs:
                self._fiber_bundles.append([])
                fb = fbs[fb]
                for f in fb:
                    f = fb[f]
                    self._fiber_bundles[-1].append(
                        np.concatenate(
                            (f['points'][:].astype(np.float32),
                             np.atleast_2d(f['radii'][:].astype(np.float32)).T),
                            axis=1))

    def ReadFiberCellFile(self, filename):
        self._fiber_bundles = []
        with h5py.File(filename, 'r') as h5f:

            fb = h5f['fiber_bundles']
            self._fiber_bundles = []
            self._fiber_bundles.append([])
            for f in fb:
                f = fb[f]
                self._fiber_bundles[-1].append(
                    Fiber(f['points'][:].flatten(), f['radii'][:]))

            self._cells_populations = []

            cells = h5f['cells/astrocytes']
            self._cells_populations.append([])
            for c in cells:
                c = cells[c]
                self._cells_populations[-1].append(
                    Cell(c['points'][:].flatten(), c['radii'][:]))

            cells = h5f['cells/oligodendrocytes']
            self._cells_populations.append([])
            for c in cells:
                c = cells[c]
                self._cells_populations[-1].append(
                    Cell(c['points'][:].flatten(), c['radii'][:]))

    def TranslateVolume(self, offset):
        if not isinstance(offset, (list, tuple)):
            raise TypeError("offset must be a list")

        for fb in self._fiber_bundles:
            for f in fb:
                f.translate(offset)

    def RotateVolume(self, phi, theta, psi):
        rot_mat = rotation.zyz(phi, theta, psi)
        for fb in self._fiber_bundles:
            for f in fb:
                f.rotate(list(rot_mat))

    def RotateVolumeAroundPoint(self, phi, theta, psi, offset):
        offset = np.array(offset)
        if offset.shape != (3,):
            raise TypeError("offset must a point")
        rot_mat = rotation.zyz(phi, theta, psi)
        for fb in self._fiber_bundles:
            for f in fb:
                f.rotate_around_point((rot_mat), offset)

    def _CheckDataLength(self):
        if (self._fiber_bundles):
            if len(self._fiber_bundles) != len(self._fiber_bundles_properties):
                raise TypeError(
                    "properties must have the same size as fiber_bundles")

        if (self._cells_populations):
            if len(self._cells_populations) != len(
                    self._cells_populations_properties):
                raise TypeError(
                    "properties must have the same size as cell_populations")

    def GenerateTissue(self, only_label=False, progress_bar=False):

        if self._dim is None:
            raise ValueError('dim not set')

        if self._dim_origin is None:
            raise ValueError('dim_origin not set')

        if self._voxel_size is None:
            raise ValueError('voxel_size not set')

        self._gen.set_volume(self._dim, self._dim_origin, self._voxel_size)
        self._CheckDataLength()
        if self._fiber_bundles:
            self._gen.set_fiber_bundles(self._fiber_bundles,
                                        self._fiber_bundles_properties)
        if self._cells_populations:
            self._gen.set_cell_populations(self._cells_populations,
                                           self._cells_populations_properties)
        label_field, vec_field, tissue_properties = self._gen.run_generation(
            only_label, progress_bar)

        return label_field, vec_field, tissue_properties

    def InitSimulation(self):
        if self._light_intensity is None:
            raise ValueError('light_intensity not set')

        if self._voxel_size is None:
            raise ValueError('voxel_size not set')

        if self._wavelength is None:
            raise ValueError('wavelength not set')

        if self._untilt_sensor is None:
            raise ValueError('untilt_sensor not set')

        if self._filter_rotations is None:
            raise ValueError('filter_rotations not set')

        self._sim.set_pli_setup(
            self._light_intensity,
            self._voxel_size,
            self._wavelength,
            self._filter_rotations,
            self._untilt_sensor,
        )

    def RunSimulation(self,
                      label_field,
                      vec_field,
                      tissue_properties,
                      theta,
                      phi,
                      step_size=1.0,
                      do_untilt=True):

        label_field = np.array(label_field, dtype=np.int32, copy=False)
        vec_field = np.array(vec_field, dtype=np.float32, copy=False)

        self.InitSimulation()

        tissue_properties = np.array(tissue_properties)
        if len(tissue_properties.shape) != 2:
            raise ValueError("tissue_properties wrong shape")

        if tissue_properties.shape[1] != 2:
            raise ValueError("tissue_properties wrong shape")

        if tissue_properties.shape[0] <= np.max(label_field.flatten()):
            raise ValueError(
                "tissue_properties.shape[0] < np.max(label_field.flatten())")

        image = self._sim.run_simulation(self._dim, label_field, vec_field,
                                         tissue_properties, theta, phi,
                                         step_size, do_untilt)
        return image

    def omp_threads(self, num_threads=2):
        ''' num_threads = 2 seems to be optimal
        '''
        if not isinstance(num_threads, int):
            raise TypeError("num_threads has to be a positiv integer")

        num_threads_gen = self._gen.set_omp_num_threads(num_threads)
        num_threads_sim = self._sim.set_omp_num_threads(num_threads)

        if num_threads_gen != num_threads_sim:
            raise ValueError("num_threads_gen != num_threads_sim")
        return num_threads_gen

    def MemoryUseage(self, unit='MB', item='all'):
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

        if item == 'label_field':
            return np.prod(self._dim) * (32 + 32) / 8 / div
        elif item == 'all':
            return np.prod(self._dim) * (32 + 32 + 3 * 32) / 8 / div
        else:
            raise ValueError('allowed is only label_field or all')

    def DimData(self):
        dim_local = self._gen.dim_local()
        dim_offset = self._gen.dim_offset()
        return dim_local, dim_offset

    def SaveAsH5(self, h5f, data, data_name, lock_dim=None):

        dim_local, dim_offset = self.DimData()

        if not isinstance(data, np.ndarray):
            raise TypeError("only numpy arrays are compatible with SaveAsH5")

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
            raise TypeError("no compatible SaveAsH5: " + data_name)

    def apply_optic(
            self,
            image_stack,
            delta_sigma=0.71,  # only for LAP!
            gain=3,  # only for LAP!
            cropping=0,  # num pixel 
            resample_mode=Image.BILINEAR):

        if self._resolution is None:
            raise TypeError("resolution is not set")

        res_image_stack = optic.apply_stack(image_stack, self._voxel_size,
                                            self._resolution, delta_sigma, gain,
                                            cropping, resample_mode)

        return res_image_stack

    def apply_resize_mask(self,
                          mask,
                          delta_sigma=0.71,
                          cropping=0,
                          resample_mode=Image.BILINEAR):
        ''' return value of mask is float, threshold has to be applyed by user
        '''

        if self._resolution is None:
            raise TypeError("resolution is not set")

        res_mask = optic.apply(np.array(mask, float), self._voxel_size,
                               self._resolution, delta_sigma, 0, cropping,
                               resample_mode)

        return res_mask

    def apply_epa(self, image_stack, mask=None):

        transmittance, direction, retardation = epa.map(image_stack)
        if mask is not None:
            transmittance[~mask] = 0
            direction[~mask] = 0
            retardation[~mask] = 0

        return transmittance, direction, retardation

    def apply_rofl(
            self,
            image_stack,
            tilt_angle=np.deg2rad(5.5),  # only LAP!
            gain=3,  # only LAP!
            mask=None,
            num_threads=2):

        rofl_direction, rofl_incl, rofl_t_rel, dirdevmap, incldevmap, treldevmap, funcmap, itermap = rofl.map(
            image_stack, tilt_angle, gain, mask, num_threads)

        return rofl_direction, rofl_incl, rofl_t_rel, (dirdevmap, incldevmap,
                                                       treldevmap, funcmap,
                                                       itermap)
