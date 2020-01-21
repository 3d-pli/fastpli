from .__generation import _Generator
from .__simulation import _Simulator

from ..version import __version__
from ..analysis import rofl, epa, affine_transformation
from ..analysis.images import fom_hsv_black
from ..simulation import optic
from ..tools import rotation
from ..tools.helper import pip_freeze

import numpy as np
import warnings

# import importlib
# if importlib.util.find_spec("tqdm"):
#     from tqdm import tqdm
# else:
#     tqdm = lambda i: i


class Simpli:
    __isfrozen = False

    def __setattr__(self, key, value):
        if self.__isfrozen and not hasattr(self, key):
            raise TypeError("%r is a frozen class" % self)
        object.__setattr__(self, key, value)

    def __freeze(self):
        self.__isfrozen = True

    def __init__(self, mpi_comm=None):

        # LIBRARIES
        self.__gen = _Generator()
        self.__sim = _Simulator()
        if mpi_comm:
            from mpi4py import MPI
            self.__gen.set_mpi_comm(MPI._addressof(mpi_comm))
            self.__sim.set_mpi_comm(MPI._addressof(mpi_comm))

        # DIMENSIONS
        self._dim = None
        self._dim_origin = np.array([0, 0, 0], dtype=float)
        self._voxel_size = None

        # GENERATION
        self._cells_populations = None
        self._cells_populations_properties = None
        self._fiber_bundles = None
        self._fiber_bundles_properties = None

        # SIMULATION
        self._filter_rotations = None
        self._flip_z_beam = False
        self._interpolate = True
        self._light_intensity = None
        self._resolution = None
        self._step_size = 1.0
        self._tilts = [(0, 0), (5, np.deg2rad(0)), (5, np.deg2rad(90)),
                       (5, np.deg2rad(180)), (5, np.deg2rad(270))]
        self._tissue_refrection = 1
        self._untilt_sensor_view = True
        self._wavelength = None

        # ANALYSIS
        self._sensor_gain = None

        # OTHER
        self._omp_num_threads = 1
        self._debug = False

        # freeze class
        self.__freeze()

    def print(self, msg):
        if self._verbose:
            print(msg)

    def get_dict(self):
        '''
        Get all member variables which are properties
        '''

        members = dict()
        for key, value in self.__dict__.items():
            if key == '_cells_populations' or key == '_fiber_bundles':
                continue
            if key.startswith("_") and not key.startswith(
                    "__") and not key.startswith("_Simpli"):
                members[key[1:]] = value
        return members

    def set_dict(self, input):
        for key, value in input.items():
            if key.startswith("_"):
                raise ValueError("member variable cant be set directly")

            if value is not None:
                setattr(self, key, value)
            else:
                warnings.warn("None value in dict detected.")

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
        dim = np.array(dim)
        if dim.dtype != int:
            raise TypeError("dim is not np.array(int)")

        if dim.size != 3:
            raise TypeError("dim.size != 3")

        if dim.ndim != 1:
            raise TypeError("dim.ndim != 1")

        self._dim = dim

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
        if voxel_size <= 0:
            raise ValueError("voxel_size <= 0")
        self._voxel_size = float(voxel_size)

    @property
    def resolution(self):
        return (self._resolution)

    @resolution.setter
    def resolution(self, resolution):
        if resolution <= 0:
            raise ValueError("resolution <= 0")
        self._resolution = float(resolution)

    def get_voi(self):
        if self._voxel_size is None:
            raise ValueError("voxel_size is not set, voi can't be calculated")

        if self._dim is None:
            raise ValueError("dim is not set, voi can't be calculated")

        if self._dim_origin is None:
            raise ValueError("dim_origin is not set, voi can't be calculated")

        voi = np.zeros((6,))
        voi[::2] = self._dim_origin
        voi[1::2] = voi[::2] + self._dim * self._voxel_size

        min = np.array(voi[0::2])
        max = np.array(voi[1::2])
        return min, max

    def set_voi(self, min, max):
        '''
        min: [x_min, y_min, z_min]
        max: [x_max, y_max, z_max]
        '''

        min = np.array(min, dtype=float)
        max = np.array(max, dtype=float)

        if min.ndim != 1 or max.ndim != 1:
            raise TypeError("min,max : ndim != 1")
        if min.size != 3 or max.size != 3:
            raise TypeError("min,max : size != 1")
        if np.all(min >= max):
            raise ValueError("min >= max")

        if self._voxel_size is None:
            raise TypeError("voxel_size is not set yet")

        self._dim = np.array((max - min) / self._voxel_size, dtype=int)
        self._dim_origin = min

    @property
    def filter_rotations(self):
        return (self._filter_rotations)

    @filter_rotations.setter
    def filter_rotations(self, filter_rotations):
        filter_rotations = np.array(filter_rotations, dtype=float)

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
    def tilts(self):
        return self._tilts

    @tilts.setter
    def tilts(self, tilts):
        ''' [[theta, phi], ...]
        '''

        if not isinstance(tilts, (list, tuple, np.ndarray)):
            raise TypeError("tilts is not a list or array")

        tilts = np.array(tilts)
        if tilts.ndim != 2:
            raise TypeError("tilts.shape != nx2")
        if tilts.shape[-1] != 2:
            raise TypeError("tilts.shape != nx2")

        self._tilts = tilts

    @property
    def sensor_gain(self):
        return self._sensor_gain

    @sensor_gain.setter
    def sensor_gain(self, sensor_gain):

        if not isinstance(sensor_gain, (int, float)):
            raise TypeError("sensor_gain is not a number")

        if sensor_gain <= 0:
            raise ValueError("sensor_gain is <= 0")

        self._sensor_gain = sensor_gain

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

        if not fbs:
            self._fiber_bundles = None
            return

        if not isinstance(fbs, (list, tuple)):
            raise TypeError("fbs is not a list")

        for fb_i, fb in enumerate(fbs):
            if not isinstance(fb, (list, tuple)):
                raise TypeError("fb is not a list")

            for f_i, f in enumerate(fb):
                fbs[fb_i][f_i] = np.array(f, dtype=float)

                if fbs[fb_i][f_i].ndim != 2:
                    raise TypeError("fiber.shape != nx4")

                if fbs[fb_i][f_i].shape[1] != 4:
                    raise TypeError("fiber.shape != nx4")

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
            raise TypeError("properties != list(list(tuples))")

        self._fiber_bundles_properties = []
        for prop in bundle_layer_properties:
            if not isinstance(prop, (list, tuple)):
                raise TypeError("properties != list(list(tuples))")

            self._fiber_bundles_properties.append([])
            for ly in prop:
                if len(ly) != 4:
                    raise TypeError("layer != (float, float, float, char)")

                if ly[1] < 0 and ly[-1] is 'r':
                    warnings.warn("birefringence negative and radial")
                if ly[1] > 0 and ly[-1] is 'p':
                    warnings.warn("birefringence positive and parallel")

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
            raise TypeError("cells_populations != list")

        for cp_i, cp in enumerate(cps):
            if not isinstance(cp, (list, tuple)):
                raise TypeError("cells_population != list")

            for c_i, c in enumerate(cp):
                cps[cp_i][c_i] = np.array(c, dtype=float)

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
                raise TypeError("len(fiber_bundles) != len(properties)")

        if self._cells_populations:
            if len(self._cells_populations) != len(
                    self._cells_populations_properties):
                raise TypeError("len(cell_populations) != len(cell_properties)")

    def generate_tissue(self, only_label=False, progress_bar=False):

        self._check_volume_input()
        self._check_generation_input()

        self.__gen.set_volume(self._dim, self._dim_origin, self._voxel_size)
        if self._fiber_bundles:
            self.__gen.set_fiber_bundles(self._fiber_bundles,
                                         self._fiber_bundles_properties)
        if self._cells_populations:
            self.__gen.set_cell_populations(self._cells_populations,
                                            self._cells_populations_properties)
        label_field, vec_field, tissue_properties = self.__gen.run_generation(
            only_label, progress_bar)

        return label_field, vec_field, tissue_properties

    def _init_pli_setup(self):
        self._check_simulation_input()
        self.__sim.set_pli_setup(self._step_size, self._light_intensity,
                                 self._voxel_size, self._wavelength,
                                 self._tissue_refrection, self._interpolate,
                                 self._untilt_sensor_view, self._flip_z_beam,
                                 self._filter_rotations)

    def _check_volume_input(self):
        if self._dim is None:
            raise ValueError("dim not set")

        if self._dim_origin is None:
            raise ValueError("dim_origin not set")

        if self._voxel_size is None:
            raise ValueError("voxel_size not set")

        self.get_voi()

    def _check_generation_input(self):
        if not self._fiber_bundles and not self._cells_populations:
            raise ValueError("fiber_bundles and cells_populations are not set")
        self._check_property_length()

    def _check_simulation_input(self):
        if self._step_size <= 0:
            raise ValueError("step_size <= 0")

        if self._light_intensity is None:
            raise ValueError("light_intensity not set")

        if self._voxel_size is None:
            raise ValueError("voxel_size not set")

        if self._wavelength is None:
            raise ValueError("wavelength not set")

        if self._filter_rotations is None:
            raise ValueError("filter_rotations not set")

    def _check_analysis_input(self):
        if not np.all(self._tilts[:,-1]==np.deg2rad(np.array([0, 0, 90, 180, 270]))):
            raise ValueError("Currently only the std tilting angles can be analysed.")

        if not self._sensor_gain:
            raise ValueError("sensor_gain not set")

    def run_simulation(self, label_field, vec_field, tissue_properties, theta,
                       phi):

        self._check_volume_input()
        self._init_pli_setup()

        tissue_properties = np.array(tissue_properties)
        if tissue_properties.ndim != 2:
            raise ValueError("tissue_properties.ndim !=2")

        if tissue_properties.shape[1] != 2:
            raise ValueError("tissue_properties.shape[1] != 2")

        label_field_ = np.array(label_field, dtype=np.int32, copy=False)
        vec_field_ = np.array(vec_field, dtype=np.float32, copy=False)

        if label_field_ is not label_field:
            warnings.warn("label_field is copied", UserWarning)
        if vec_field_ is not vec_field:
            warnings.warn("vec_field is copied", UserWarning)

        images = self.__sim.run_simulation(self._dim, label_field_, vec_field_,
                                           tissue_properties, theta, phi)
        if np.min(images.flatten()) < 0:
            raise ValueError("intensity < 0 detected")

        return images

    def run_simulation_pipeline(self,
                                label_field,
                                vector_field,
                                tissue_properties,
                                h5f=None):

        if self._tilts is None:
            raise ValueError("tilts not set")

        tilting_stack = [None] * len(self._tilts)

        for t, tilt in enumerate(self._tilts):
            theta, phi = tilt[0], tilt[1]
            images = self.run_simulation(label_field, vector_field,
                                         tissue_properties, theta, phi)

            if h5f:
                h5f['simulation/data/' + str(t)] = images
                h5f['simulation/data/' + str(t)].attrs['theta'] = theta
                h5f['simulation/data/' + str(t)].attrs['phi'] = phi

            # apply optic to simulation
            new_images = self.apply_optic(images, gain=self._sensor_gain)

            if h5f:
                h5f['simulation/optic/' + str(t)] = new_images
                h5f['simulation/optic/' + str(t)].attrs['theta'] = theta
                h5f['simulation/optic/' + str(t)].attrs['phi'] = phi

            # calculate modalities
            epa = self.apply_epa(new_images)

            if h5f:
                h5f['analysis/epa/' + str(t) + '/transmittance'] = epa[0]
                h5f['analysis/epa/' + str(t) + '/direction'] = np.rad2deg(
                    epa[1])
                h5f['analysis/epa/' + str(t) + '/retardation'] = epa[2]

                h5f['analysis/epa/' + str(t) +
                    '/transmittance'].attrs['theta'] = theta
                h5f['analysis/epa/' + str(t) +
                    '/transmittance'].attrs['phi'] = phi
                h5f['analysis/epa/' + str(t) +
                    '/direction'].attrs['theta'] = theta
                h5f['analysis/epa/' + str(t) + '/direction'].attrs['phi'] = phi
                h5f['analysis/epa/' + str(t) +
                    '/retardation'].attrs['theta'] = theta
                h5f['analysis/epa/' + str(t) +
                    '/retardation'].attrs['phi'] = phi

            tilting_stack[t] = new_images

        # pseudo mask
        mask = np.sum(label_field, 2) > 0
        mask = self.apply_resize(1.0 * mask) > 0.1
        h5f['simulation/optic/mask'] = np.uint8(mask)

        tilting_stack = np.array(tilting_stack)
        while tilting_stack.ndim < 4:
            tilting_stack = np.expand_dims(tilting_stack, axis=-2)

        rofl_direction, rofl_incl, rofl_t_rel, _ = self.apply_rofl(
            tilting_stack,
            tilt_angle=np.deg2rad(self._tilts[-1][0]),
            gain=self._sensor_gain,
            mask=None)

        if h5f:
            h5f['analysis/rofl/direction'] = rofl_direction
            h5f['analysis/rofl/inclination'] = rofl_incl
            h5f['analysis/rofl/trel'] = rofl_t_rel

            h5f['/parameter/dict'] = str(self.get_dict())

        fom = fom_hsv_black(rofl_direction, rofl_incl)

        return tilting_stack, (rofl_direction, rofl_incl, rofl_t_rel), fom

    def run_pipeline(self, h5f=None, script=None, save=[], verbose=False):

        self._check_volume_input()
        self._check_generation_input()
        self._check_simulation_input()
        if self._tilts is None:
            raise ValueError("tilts not set")

        if h5f:
            h5f['/parameter/version'] = __version__
            h5f['/parameter/pip_freeze'] = pip_freeze()
            if script:
                h5f['/parameter/script'] = script

        # run tissue generation
        label_field, vec_field, tissue_properties = self.generate_tissue()

        if h5f and save:

            if not "label_field" in save or not "vec_field" in save:
                raise ValueError("only label_field and vec_field are allowed")

            if "label_field" in save:
                dset = h5f.create_dataset('tissue/label_field',
                                          label_field.shape,
                                          dtype=np.uint16,
                                          compression='gzip',
                                          compression_opts=1)
                dset[:] = label_field

            if "vec_field" in save:
                dset = h5f.create_dataset('tissue/vec_field',
                                          vec_field.shape,
                                          dtype=np.float32,
                                          compression='gzip',
                                          compression_opts=1)
                dset[:] = vec_field
            


            h5f['tissue/tissue_properties'] = tissue_properties

        tilting_stack, (rofl_direction, rofl_incl,
                        rofl_t_rel), fom = self.run_simulation_pipeline(
                            label_field, vec_field, tissue_properties, h5f)

        return (label_field, vec_field,
                tissue_properties), tilting_stack, (rofl_direction, rofl_incl,
                                                    rofl_t_rel), fom

    @property
    def omp_num_threads(self):
        return self._omp_num_threads

    @omp_num_threads.setter
    def omp_num_threads(self, num_threads):

        if not isinstance(num_threads, int):
            raise TypeError("num_threads != int")

        if num_threads <= 0:
            raise TypeError("num_threads <= 0")

        num_threads_gen = self.__gen.set_omp_num_threads(num_threads)
        num_threads_sim = self.__sim.set_omp_num_threads(num_threads)

        if num_threads_gen != num_threads_sim:
            raise AssertionError("num_threads_gen != num_threads_sim")

        if num_threads_gen != num_threads:
            warnings.warn("reduced num_threads: " + str(num_threads_gen),
                          UserWarning)

        self._omp_num_threads = num_threads_gen

    def memory_usage(self, unit='MB', item='all'):

        if unit == 'MB':
            div = 1024**2
        elif unit == 'GB':
            div = 1024**3
        else:
            raise ValueError("allowed is only \"MB\", \"GB\"")

        if self._dim is None:
            raise TypeError("dimension not set yet")

        if item == 'label_field':
            # label_field + distance_array
            return np.prod(self._dim) * (32 + 32) / 8 / div
        elif item == 'all':
            return np.prod(self._dim) * (32 + 32 + 3 * 32) / 8 / div
        else:
            raise ValueError("allowed is only \"label_field\" or \"all\"")

    def save_mpi_array_as_h5(self, h5f, input, data_name, lock_dim=None):
        '''
        simpli can be seperated into different mpi processes.
        This function provides a parallel hdf5 io to save data
        inside the same h5-file.
        '''
        # TODO: check functionality

        dim_local = self.__gen.dim_local()
        dim_offset = self.__gen.dim_offset()

        if not isinstance(input, np.ndarray):
            raise TypeError(
                "only numpy arrays are compatible with save_mpi_array_as_h5")

        dset_dim = np.copy(self._dim)
        if len(input.shape) < len(dset_dim):
            dset_dim = dset_dim[:len(input.shape)]
        if len(input.shape) > len(dset_dim):
            dset_dim = np.append(dset_dim, input.shape[3:])

        if lock_dim:
            if isinstance(lock_dim, int):
                lock_dim = [lock_dim]

            lock_dim = list(lock_dim)
            for l in lock_dim:
                dset_dim[l] = input.shape[l]

        dset = h5f.create_dataset(data_name, dset_dim, dtype=input.dtype)

        if len(input.shape) == 2:
            if input.size * input.itemsize > 2 * (2**10)**3:  # 2 GB
                for i in range(input.shape[0]):
                    dset[i + dim_offset[0], dim_offset[1]:dim_offset[1] +
                         dim_local[1]] = input[i, :]
            else:
                dset[dim_offset[0]:dim_offset[0] + dim_local[0],
                     dim_offset[1]:dim_offset[1] + dim_local[1]] = input

        elif len(input.shape) == 3:
            if input.size * input.itemsize > 2 * (2**10)**3:  # 2 GB
                for i in range(input.shape[0]):
                    dset[i + dim_offset[0],
                         dim_offset[1]:dim_offset[1] + dim_local[1],
                         dim_offset[2]:dim_offset[2] +
                         dim_local[2]] = input[i, :]
            else:
                dset[dim_offset[0]:dim_offset[0] + dim_local[0],
                     dim_offset[1]:dim_offset[1] + dim_local[1],
                     dim_offset[2]:dim_offset[2] + dim_local[2]] = input

        elif len(input.shape) > 3:
            if input.size * input.itemsize > 2 * (2**10)**3:  # 2 GB
                for i in range(input.shape[0]):
                    dset[i + dim_offset[0],
                         dim_offset[1]:dim_offset[1] + dim_local[1],
                         dim_offset[2]:dim_offset[2] +
                         dim_local[2], :] = input[i, :]
            else:
                dset[dim_offset[0]:dim_offset[0] + dim_local[0],
                     dim_offset[1]:dim_offset[1] + dim_local[1],
                     dim_offset[2]:dim_offset[2] + dim_local[2], :] = input

        else:
            raise TypeError("no compatible save_mpi_array_as_h5: " + data_name)

    def apply_optic(
        self,
        input,
        delta_sigma=0.71,  # only for LAP!
        gain=3.0,  # only for LAP!
        order=1,
        resample=False,
        mp_pool=None):
        ''' 
        input = np.array(x,y(,rho))
        use resample, resize will be removed in the future
        '''

        if self._resolution is None:
            raise TypeError("resolution is not set")

        if self._voxel_size is None:
            raise TypeError("voxel_size is not set")

        input = np.atleast_3d(np.array(input))
        if input.ndim > 3:
            raise TypeError("input can be 1d, 2d or 3d")

        scale = self._voxel_size / self._resolution
        size = np.array(np.round(np.array(input.shape[0:2]) * scale), dtype=int)

        output = np.empty((size[0], size[1], input.shape[2]), dtype=input.dtype)

        if resample:
            filter_resize = lambda input, delta_sigma, scale, order: optic.filter_resample(
                input, delta_sigma.scale)
        else:
            filter_resize = optic.filter_resize

        if mp_pool:
            chunk = [(input[:, :, i], delta_sigma, scale, order)
                     for i in range(input.shape[2])]
            results = mp_pool.starmap(filter_resize, chunk)
            for i in range(input.shape[2]):
                output[:, :, i] = results[i]
        else:
            for i in range(input.shape[2]):
                output[:, :, i] = filter_resize(input[:, :, i], delta_sigma,
                                                scale, order)

        if np.min(output.flatten()) < 0:
            raise AssertionError("intensity < 0 detected")

        if gain > 0:
            output = optic.add_noise(output, gain)

        return np.squeeze(output)

    def apply_resize(self, input, order=1, mp_pool=None):
        ''' 
        input = np.array(x,y(,rho))
        '''

        if self._resolution is None:
            raise TypeError("resolution is not set")

        if self._voxel_size is None:
            raise TypeError("voxel_size is not set")

        input = np.atleast_3d(np.array(input))
        if input.ndim > 3:
            raise TypeError("input can be 1d, 2d or 3d")

        scale = self._voxel_size / self._resolution
        size = np.array(np.round(np.array(input.shape[0:2]) * scale), dtype=int)

        output = np.empty((size[0], size[1], input.shape[2]), dtype=input.dtype)
        if mp_pool:
            input_data = [
                (input[:, :, i], scale, order) for i in range(input.shape[2])
            ]
            results = mp_pool.starmap(optic.resize, input_data)
            for i in range(input.shape[2]):
                output[:, :, i] = results[i]
        else:
            for i in range(input.shape[2]):
                output[:, :, i] = optic.resize(input[:, :, i], scale, order)
        return np.squeeze(output)

    def apply_untilt(self, input, theta, phi, mode='nearest'):
        ''' 
        input = np.array(x,y(,rho))
        '''

        if theta == 0:
            return input

        # calculate transformation matrix
        p = self._dim.copy()
        p_rot = 0.5 * np.array([p[0], p[1], 0])
        p_out = np.array([[p[0], p[1], 0], [p[0], 0, 0], [0, p[1], 0]])
        rot = rotation.theta_phi(-theta, phi)

        p_in = np.array([np.dot(rot, p - p_rot) + p_rot for p in p_out])

        # TODO: refraction has to be implemented
        M = affine_transformation.calc_matrix(p_in[:, :2], p_out[:, :2])

        output = affine_transformation.image(input, M, mode)
        if input.ndim == 3:
            output = np.atleast_3d(output)

        return output

    def apply_epa(self, input, mask=None):
        ''' 
        input = np.array(x,y,rho)
        '''

        transmittance, direction, retardation = epa.epa(input)
        if mask is not None:
            transmittance[np.invert(mask)] = float('nan')
            direction[np.invert(mask)] = float('nan')
            retardation[np.invert(mask)] = float('nan')

        return transmittance, direction, retardation

    def apply_rofl(
        self,
        input,
        tilt_angle=np.deg2rad(5.5),  # only LAP!
        gain=3.0,  # only LAP!
        dir_offset=0,
        mask=None,
        mp_pool=None,
        grad_mode=False):
        ''' 
        input = np.array(tilt,x,y,rho)
        '''

        input = np.array(input, copy=False)

        if input.ndim != 4:
            raise TypeError("input: np.array([tilts,x,y,stack])")

        if input.shape[0] != 5:
            raise ValueError("input need 1 + 4 measurements")

        if input.shape[-1] <= 3:
            raise ValueError("input needs at least 3 equidistand rotations")

        if gain <= 0:
            raise ValueError("rofl gain <= 0")

        if mask is None:
            mask = np.ones((input.shape[1], input.shape[2]), bool)

        directionmap = np.empty_like(mask, dtype=input.dtype)
        inclmap = np.empty_like(mask, dtype=input.dtype)
        trelmap = np.empty_like(mask, dtype=input.dtype)
        dirdevmap = np.empty_like(mask, dtype=input.dtype)
        incldevmap = np.empty_like(mask, dtype=input.dtype)
        treldevmap = np.empty_like(mask, dtype=input.dtype)
        funcmap = np.empty_like(mask, dtype=input.dtype)
        itermap = np.empty_like(mask, dtype=input.dtype)

        if mp_pool:
            for j in range(input.shape[2]):
                chunk = [(input[:, i,
                                j, :], tilt_angle, gain, dir_offset, grad_mode)
                         for i in range(input.shape[1])]
                results = mp_pool.starmap(rofl.rofl, chunk)

                for i, result in enumerate(results):
                    directionmap[i, j], inclmap[i, j], trelmap[i, j], dirdevmap[
                        i, j], incldevmap[i, j], treldevmap[i, j], funcmap[
                            i, j], itermap[i, j] = result
        else:
            for i in range(input.shape[1]):
                for j in range(input.shape[2]):
                    directionmap[i, j], inclmap[i, j], trelmap[i, j], dirdevmap[
                        i, j], incldevmap[i, j], treldevmap[i, j], funcmap[
                            i,
                            j], itermap[i,
                                        j] = rofl.rofl(input[:, i, j, :],
                                                       tilt_angle, gain,
                                                       dir_offset, grad_mode)

        return directionmap, inclmap, trelmap, (dirdevmap, incldevmap,
                                                treldevmap, funcmap, itermap)
