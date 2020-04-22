# -*- coding: utf-8 -*-
"""
Simpli Class
"""

from .__generation import _Generator
from .__simulation import _Simulator

from . import optic
from .. import analysis
from .. import objects
from .. import tools
from .. import version

import numpy as np
import warnings


class Simpli:
    """
    Simpli Class for simulating 3D-PLI images
    """

    __is_frozen = False

    def __setattr__(self, key, value):
        if self.__is_frozen and not hasattr(self, key):
            raise TypeError(f"{self.__class__.__name__} is a frozen class")
        object.__setattr__(self, key, value)

    def __freeze(self):
        self.__is_frozen = True

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
        self._pixel_size = None
        self._optical_sigma = None
        self._step_size = 1.0
        self._tilts = None
        self._tissue_refrection = 1
        self._untilt_sensor_view = True
        self._wavelength = None

        # ANALYSIS
        self._sensor_gain = None

        # OTHER
        self._omp_num_threads = 1
        self._verbose = False

        # freeze class
        self.__freeze()

    def _print(self, msg):
        if self._verbose:
            print(msg)

    def get_dict(self):
        """ Get all member variables which are properties """
        members = dict()
        for key, value in self.__dict__.items():
            if key == '_cells_populations' or key == '_fiber_bundles':
                continue
            if key.startswith("_") and not key.startswith(
                    "__") and not key.startswith("_Simpli"):

                if isinstance(value, np.ndarray):
                    members[key[1:]] = value.tolist()
                else:
                    members[key[1:]] = value

        return members

    def set_dict(self, input):
        """ Set dictionary of variables to class members """
        for key, value in input.items():
            if key.startswith("_"):
                raise ValueError("member variable cant be set directly")

            if value is not None:
                setattr(self, key, value)
            else:
                warnings.warn("None value in dict detected")

    @property
    def dim(self):
        """ dimension of volume in voxel: (3)-array """
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
        """ origin of volume in µm: (3)-array """
        return self._dim_origin

    @dim_origin.setter
    def dim_origin(self, dim_origin):
        self._dim_origin = np.array(dim_origin, dtype=float)

    @property
    def voxel_size(self):
        """ size of voxel in µm: float """
        return (self._voxel_size)

    @voxel_size.setter
    def voxel_size(self, voxel_size):
        if voxel_size <= 0:
            raise ValueError("voxel_size <= 0")

        flag = False
        if self._voxel_size is not None and self._dim is not None and self._dim_origin is not None:
            min, max = self.get_voi()
            flag = True

        self._voxel_size = float(voxel_size)

        if flag and self._dim is not None and self._dim_origin is not None:
            self.set_voi(min, max)
            self._print("voxel_size: recalculated dimensions")

    @property
    def pixel_size(self):
        """ pixel size of resulting optical image in µm: float """
        return (self._pixel_size)

    @pixel_size.setter
    def pixel_size(self, pixel_size):
        if pixel_size <= 0:
            raise ValueError("pixel_size <= 0")
        self._pixel_size = float(pixel_size)

    def get_voi(self):
        """ get volume of interest """
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
        """ 
        set volume of interest

        min: [x_min, y_min, z_min]
        max: [x_max, y_max, z_max]
        """

        min = np.array(min, dtype=float)
        max = np.array(max, dtype=float)

        if min.ndim != 1 or max.ndim != 1:
            raise TypeError("min,max : ndim != 1")
        if min.size != 3 or max.size != 3:
            raise TypeError("min,max : size != 1")
        if np.any(min >= max):
            raise ValueError("min >= max")

        if self._voxel_size is None:
            raise TypeError("voxel_size is not set yet")

        self._dim = np.array(np.round((max - min) / self._voxel_size),
                             dtype=int)
        self._dim_origin = min

    @property
    def filter_rotations(self):
        """ list of filter rotation position in radiant: [float] """
        return (self._filter_rotations)

    @filter_rotations.setter
    def filter_rotations(self, filter_rotations):
        filter_rotations = np.array(filter_rotations, dtype=float)

        if filter_rotations.size == 0 or filter_rotations.ndim != 1:
            raise TypeError("filter_rotations : nx1")

        self._filter_rotations = filter_rotations

    @property
    def light_intensity(self):
        """ initial light intensity value in a.u.: float """
        return self._light_intensity

    @light_intensity.setter
    def light_intensity(self, light_intensity):
        self._light_intensity = float(light_intensity)

    @property
    def tilts(self):
        """ list of spherical tilting angles [[theta, phi], ...] in radiant:
        [[float, float], ...] 
        """
        return self._tilts

    @tilts.setter
    def tilts(self, tilts):
        if not isinstance(tilts, (list, tuple, np.ndarray)):
            raise TypeError("tilts is not a list or array")

        tilts = np.array(tilts, ndmin=2)
        if tilts.ndim != 2:
            raise TypeError("tilts.shape != nx2")
        if tilts.shape[-1] != 2:
            raise TypeError("tilts.shape != nx2")

        self._tilts = tilts

    @property
    def sensor_gain(self):
        """ sensor gain value for calculating ccd camera noise in a.u.: float """
        return self._sensor_gain

    @sensor_gain.setter
    def sensor_gain(self, sensor_gain):

        if not isinstance(sensor_gain, (int, float)):
            raise TypeError("sensor_gain is not a number")

        if sensor_gain <= 0:
            raise ValueError("sensor_gain is <= 0")

        self._sensor_gain = sensor_gain

    @property
    def optical_sigma(self):
        """ optical sigma for applying a gaussian convolution to the image
        after resizing in pixel_size: float 
        """
        return self._optical_sigma

    @optical_sigma.setter
    def optical_sigma(self, optical_sigma):

        if not isinstance(optical_sigma, (int, float)):
            raise TypeError("optical_sigma is not a number")

        if optical_sigma <= 0:
            raise ValueError("optical_sigma is <= 0")

        self._optical_sigma = optical_sigma

    @property
    def wavelength(self):
        """ wavelength of light in nm: float """
        return self._wavelength

    @wavelength.setter
    def wavelength(self, wavelength):
        self._wavelength = float(wavelength)

    @property
    def tissue_refrection(self):
        """ tissue refractive index in a.u.: float """
        return self._tissue_refrection

    @tissue_refrection.setter
    def tissue_refrection(self, tissue_refrection):
        self._tissue_refrection = float(tissue_refrection)

    @property
    def step_size(self):
        """ step size for light in voxel: float """
        return self._step_size

    @step_size.setter
    def step_size(self, step_size):
        self._step_size = float(step_size)

    @property
    def interpolate(self):
        """ (de)activate interpolation of vectors in simulation: bool """
        return self._interpolate

    @interpolate.setter
    def interpolate(self, interpolate):
        self._interpolate = bool(interpolate)

    @property
    def verbose(self):
        """ additional information will be printed: bool """
        return self._verbose

    @verbose.setter
    def verbose(self, verbose):
        self._verbose = bool(verbose)

    @property
    def untilt_sensor_view(self):
        """ untilt the image by adapted initial light coordinates: bool
        
        otherwise the image has to be untilted with an affine transformation
        """
        return self._untilt_sensor_view

    @untilt_sensor_view.setter
    def untilt_sensor_view(self, untilt_sensor_view):
        self._untilt_sensor_view = bool(untilt_sensor_view)

    @property
    def flip_z_beam(self):
        """ flip the direction of light along the z-axis e.g. PM: bool """
        return self._flip_z_beam

    @flip_z_beam.setter
    def flip_z_beam(self, flip_z_beam):
        self._flip_z_beam = bool(flip_z_beam)

    @property
    def fiber_bundles(self):
        """ set fiber_bundles: [[(,4)-array]] """
        return self._fiber_bundles

    @fiber_bundles.setter
    def fiber_bundles(self, fbs):
        fbs = objects.fiber_bundles.Cast(fbs)
        self._fiber_bundles = fbs

    @property
    def fiber_bundles_properties(self):
        """ set physical properties of fiber_bundles per fiber_bundle [(dn, mu)]:
        [(float, float)]
        """
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
        """ set cells_populations: [[(,4)-array]] """
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

                if cps[cp_i][c_i].ndim != 2:
                    raise TypeError("cell size need to be nx4")

                if cps[cp_i][c_i].shape[1] != 4:
                    raise TypeError("cell size need to be nx4")

        self._cells_populations = cps

    @property
    def cells_populations_properties(self):
        """ set physical properties of fiber_bundles per fiber_bundle [(dn, mu)]:
        [(float, float)]
        """
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
                raise TypeError("cell properties must be a list of 2 arguments")

            if len(prop) != 2:
                raise TypeError("cell properties must be a list of 2 arguments")

        self._cells_populations_properties = cells_populations_properties

    def _check_property_length(self):
        if self._fiber_bundles:
            if len(self._fiber_bundles) != len(self._fiber_bundles_properties):
                raise TypeError(
                    "len(fiber_bundles) != len(fiber_bundles_properties)\n\
                        For each fiber_bundle there hast to be a [(prop), ...]")

        if self._cells_populations:
            if len(self._cells_populations) != len(
                    self._cells_populations_properties):
                raise TypeError("len(cell_populations) != len(cell_properties)")

    def generate_tissue(self, only_label=False, progress_bar=False):
        """ generating discret tissue for simulation
        """

        self._print("Generate Tissue")
        self._print(f"Memory needed: ~{np.ceil(self.memory_usage()):.0f} MB")

        self._check_volume_input()
        self._check_generation_input()

        self.__gen.set_volume(self._dim, self._dim_origin, self._voxel_size)
        if self._fiber_bundles:
            self.__gen.set_fiber_bundles(self._fiber_bundles,
                                         self._fiber_bundles_properties)
        if self._cells_populations:
            self.__gen.set_cell_populations(self._cells_populations,
                                            self._cells_populations_properties)
        tissue, optical_axis, tissue_properties = self.__gen.run_generation(
            only_label, progress_bar)

        if not np.any(tissue):
            warnings.warn(
                "All labels are 0. Usually this means, that the VOI contains no\
                     fiber_bundles or that the voxel_size is to large.",
                UserWarning)

        return tissue, optical_axis, tissue_properties

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

    def run_simulation(self, tissue, optical_axis, tissue_properties, theta,
                       phi):
        """ running simulation. Input from tissue_generation

        tissue and optical_axis will be passed by reference
        theta, phi: tilting angle in radiant
        """

        self._check_volume_input()
        self._init_pli_setup()

        tissue_properties = np.array(tissue_properties)
        if tissue_properties.ndim != 2:
            raise ValueError("tissue_properties.ndim !=2")

        if tissue_properties.shape[1] != 2:
            raise ValueError("tissue_properties.shape[1] != 2")

        tissue_ = np.array(tissue, dtype=np.int32, copy=False)
        optical_axis_ = np.array(optical_axis, dtype=np.float32, copy=False)

        if tissue_ is not tissue:
            warnings.warn("tissue is copied", UserWarning)
        if optical_axis_ is not optical_axis:
            warnings.warn("optical_axis is copied", UserWarning)

        images = self.__sim.run_simulation(self._dim, tissue_, optical_axis_,
                                           tissue_properties, theta, phi)
        if np.min(images.flatten()) < 0:
            raise ValueError("intensity < 0 detected")

        return images

    def save_parameter_h5(self, h5f, script=None):
        """ Saves class members without fiber_bundles in hdf5 file. """
        h5f.attrs['fastpli/simpli'] = str(self.get_dict())
        h5f.attrs['fastpli/version'] = version.__version__
        h5f.attrs['fastpli/compiler'] = version.__compiler__
        h5f.attrs['fastpli/libraries'] = version.__libraries__
        h5f.attrs['fastpli/pip_freeze'] = tools.helper.pip_freeze()
        if script:
            h5f.attrs['script'] = script

    def run_tissue_pipeline(self,
                            h5f=None,
                            script=None,
                            save=['tissue', 'optical_axis']):
        """ Automatic pipeline for tissue generation with save options """

        self._check_volume_input()
        self._check_generation_input()

        if 'all' in save or 'simulation' in save:
            save = save + ['tissue', 'optical_axis']

        # run tissue generation
        tissue, optical_axis, tissue_properties = self.generate_tissue()

        if h5f and "tissue" in save:
            self._print("Save tissue")
            dset = h5f.create_dataset('tissue/tissue',
                                      tissue.shape,
                                      dtype=np.uint16,
                                      compression='gzip',
                                      compression_opts=1)
            dset[:] = tissue

        if h5f and "optical_axis" in save:
            self._print("Save optical_axis")
            dset = h5f.create_dataset('tissue/optical_axis',
                                      optical_axis.shape,
                                      dtype=np.float32,
                                      compression='gzip',
                                      compression_opts=1)
            dset[:] = optical_axis

        if h5f and "tissue" in save:
            h5f['tissue/properties'] = tissue_properties

        return tissue, optical_axis, tissue_properties

    def crop_tilt_pixel(self):
        """ crop affected boundary pixel from tilted images """
        if not self.untilt_sensor_view:
            raise ValueError("currently only for untilt_sensor_view=True")
        if self._voxel_size is None:
            raise ValueError("voxel_size is not set")
        if self._pixel_size is None:
            raise ValueError("pixel_size is not set")
        if self._tilts is None:
            raise ValueError("tilts are not set")

        delta_pixel = int(
            np.ceil(
                np.tan(np.max(np.abs(self._tilts[:, 0]))) * self._dim[-1] / 2 *
                self._voxel_size / self._pixel_size))
        return delta_pixel

    def crop_tilt_voxel(self):
        """ get number of affected boundary voxel from tilted images """
        delta_voxel = int(
            np.round(self.crop_tilt_pixel() * self._pixel_size /
                     self._voxel_size))
        return delta_voxel

    def add_crop_tilt_halo(self):
        """ add number of necessary boundary voxel from tilted images """
        self.dim_origin[:2] -= self.crop_tilt_voxel() * self.voxel_size
        self.dim[:2] += 2 * self.crop_tilt_voxel()

    def rm_crop_tilt_halo(self, input):
        """ remove number of added boundary voxel from tilted images """
        delta_voxel = self.crop_tilt_voxel()
        if delta_voxel == 0:
            return input
        return input[delta_voxel:-delta_voxel, delta_voxel:-delta_voxel, :]

    def run_simulation_pipeline(self,
                                tissue,
                                optical_axis,
                                tissue_properties,
                                h5f=None,
                                save=['data', 'optic', 'epa', 'mask', 'rofl'],
                                crop_tilt=False,
                                mp_pool=None):
        """ Automatic pipeline for simulation and analysis with save options """

        if 'all' in save or 'simulation' in save:
            save = save + ['data', 'optic', 'epa', 'mask', 'rofl']

        if self._tilts is None:
            raise ValueError("tilts is not set")
        if self._optical_sigma is None:
            raise ValueError("optical_sigma is not set")

        flag_rofl = True
        if np.any(self._tilts[:, 1] != np.deg2rad([0, 0, 90, 180, 270])
                 ) or self._tilts[0, 0] != 0 or np.any(
                     self._tilts[1:, 0] != self._tilts[1, 0]):
            warnings.warn("Tilts not suitable for ROFL. Skipping analysis")
            flag_rofl = False

        tilting_stack = [None] * len(self._tilts)

        self._print("Simulate tilts:")
        for t, tilt in enumerate(self._tilts):
            theta, phi = tilt[0], tilt[1]
            self._print(
                f"Tilt {t}: theta: {np.rad2deg(theta):.1f} deg, phi: {np.rad2deg(phi):.1f} deg"
            )
            images = self.run_simulation(tissue, optical_axis,
                                         tissue_properties, theta, phi)

            if crop_tilt:
                images = self.rm_crop_tilt_halo(images)

            if h5f and 'data' in save:
                self._print("Save data")
                h5f['simulation/data/' + str(t)] = images
                h5f['simulation/data/' + str(t)].attrs['theta'] = theta
                h5f['simulation/data/' + str(t)].attrs['phi'] = phi

            # apply optic to simulation
            new_images = self.apply_optic(images, mp_pool=mp_pool)

            if h5f and 'optic' in save:
                self._print("Save optic")
                h5f['simulation/optic/' + str(t)] = new_images
                h5f['simulation/optic/' + str(t)].attrs['theta'] = theta
                h5f['simulation/optic/' + str(t)].attrs['phi'] = phi

            # calculate modalities
            epa = self.apply_epa(new_images)

            if h5f and 'epa' in save:
                self._print("Save epa")
                h5f['analysis/epa/' + str(t) + '/transmittance'] = epa[0]
                h5f['analysis/epa/' + str(t) + '/direction'] = epa[1]
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
        mask = np.sum(tissue, 2) > 0
        mask = self.apply_optic_resample(1.0 * mask, mp_pool=mp_pool) > 0.1
        if h5f and 'mask' in save:
            self._print("Save mask")
            h5f['simulation/optic/mask'] = np.uint8(mask)

        tilting_stack = np.array(tilting_stack)
        while tilting_stack.ndim < 4:
            tilting_stack = np.expand_dims(tilting_stack, axis=-2)

        if flag_rofl:
            rofl_direction, rofl_incl, rofl_t_rel, (
                rofl_direction_conf, rofl_incl_conf, rofl_t_rel_conf, rofl_func,
                rofl_n_iter) = self.apply_rofl(tilting_stack,
                                               mask=None,
                                               mp_pool=mp_pool)
        else:
            rofl_direction = None
            rofl_incl = None
            rofl_t_rel = None

            rofl_direction_conf = None
            rofl_incl_conf = None
            rofl_t_rel_conf = None
            rofl_func = None
            rofl_n_iter = None

        if h5f and flag_rofl and 'rofl' in save:
            self._print("Save rofl")
            h5f['analysis/rofl/direction'] = rofl_direction
            h5f['analysis/rofl/inclination'] = rofl_incl
            h5f['analysis/rofl/t_rel'] = rofl_t_rel

            h5f['analysis/rofl/direction_conf'] = rofl_direction_conf,
            h5f['analysis/rofl/inclination_conf'] = rofl_incl_conf,
            h5f['analysis/rofl/t_rel_conf'] = rofl_t_rel_conf,
            h5f['analysis/rofl/func'] = rofl_func,
            h5f['analysis/rofl/n_iter'] = rofl_n_iter

        if flag_rofl:
            fom = analysis.images.fom_hsv_black(rofl_direction, rofl_incl)
        else:
            fom = None

        return tilting_stack, (rofl_direction, rofl_incl, rofl_t_rel), fom

    def run_pipeline(
        self,
        h5f=None,
        script=None,
        save=['tissue', 'optical_axis', 'data', 'optic', 'epa', 'mask', 'rofl'],
        crop_tilt=False,
        mp_pool=None):
        """ Automatic tissue generation and simulation pipeline with save options """

        if 'all' in save:
            save = [
                'tissue', 'optical_axis', 'data', 'optic', 'epa', 'mask', 'rofl'
            ]
        if 'tissue' in save:
            save = save + ['tissue', 'optical_axis']
        if 'simulation' in save:
            save = save + ['data', 'optic', 'epa', 'mask', 'rofl']

        # check parameters
        self._check_volume_input()
        self._check_generation_input()
        self._check_simulation_input()
        if self._tilts is None:
            raise ValueError("tilts is not set")
        if self._optical_sigma is None:
            raise ValueError("optical_sigma is not set")

        # save parameters
        if h5f:
            self.save_parameter_h5(h5f, script)

        # run tissue generation
        tissue, optical_axis, tissue_properties = self.run_tissue_pipeline(
            h5f=h5f, script=script, save=save)

        # run simulation and analysis
        tilting_stack, (rofl_direction, rofl_incl,
                        rofl_t_rel), fom = self.run_simulation_pipeline(
                            tissue,
                            optical_axis,
                            tissue_properties,
                            h5f=h5f,
                            save=save,
                            crop_tilt=crop_tilt,
                            mp_pool=mp_pool)

        return (tissue, optical_axis,
                tissue_properties), tilting_stack, (rofl_direction, rofl_incl,
                                                    rofl_t_rel), fom

    @property
    def omp_num_threads(self):
        """ get/set number of omp threads """
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
        """ print expected memory ussage """
        if unit == 'MB':
            div = 1024**2
        elif unit == 'GB':
            div = 1024**3
        else:
            raise ValueError("allowed is only \"MB\", \"GB\"")

        if self._dim is None:
            raise TypeError("dimension not set yet")

        if item == 'tissue':
            # tissue + distance_array
            return np.prod(self._dim) * (32 + 32) / 8 / div
        elif item == 'all':
            return np.prod(self._dim) * (32 + 32 + 3 * 32) / 8 / div
        else:
            raise ValueError("allowed is only \"tissue\" or \"all\"")

    def save_mpi_array_as_h5(self, h5f, input, data_name, lock_dim=None):
        """
        simpli can be seperated into different mpi processes.
        This function provides a parallel hdf5 io to save data
        inside the same h5-file.
        """
        # TODO: check functionality

        dim_local = self.__gen.dim_local()
        dim_offset = self.__gen.dim_offset()
        input = np.array(input, copy=False)

        if not isinstance(input, np.ndarray):
            raise TypeError(
                "only numpy arrays are compatible with save_mpi_array_as_h5")

        dset_dim = np.copy(self._dim)
        if input.ndim < len(dset_dim):
            dset_dim = dset_dim[:len(input.shape)]
        if input.ndim > len(dset_dim):
            dset_dim = np.append(dset_dim, input.shape[3:])

        if lock_dim:
            if isinstance(lock_dim, int):
                lock_dim = [lock_dim]

            lock_dim = list(lock_dim)
            for l in lock_dim:
                dset_dim[l] = input.shape[l]

        dset = h5f.create_dataset(data_name, dset_dim, dtype=input.dtype)

        if input.ndim == 2:
            if input.size * input.itemsize > 2 * (2**10)**3:  # 2 GB
                for i in range(input.shape[0]):
                    dset[i + dim_offset[0], dim_offset[1]:dim_offset[1] +
                         dim_local[1]] = input[i, :]
            else:
                dset[dim_offset[0]:dim_offset[0] + dim_local[0],
                     dim_offset[1]:dim_offset[1] + dim_local[1]] = input

        elif input.ndim == 3:
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

        elif input.ndim > 3:
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

    def apply_optic(self, input, shift=(0, 0), mp_pool=None):
        """ applies optical resample, convolution and noise to image
        input: np.array(x,y(,rho))
        """

        if self._sensor_gain is None:
            raise ValueError("sensor_gain not set")

        output = self.apply_optic_resample(input, mp_pool=mp_pool)

        if np.amin(output) < 0:
            raise AssertionError("intensity < 0 detected")

        if self._sensor_gain > 0:
            output = optic.add_noise(output, self._sensor_gain)

        return output

    def apply_optic_resample(self, input, shift=(0, 0), mp_pool=None):
        """ applies optical resample, convolution and noise to image
        input: np.array(x,y(,rho))
        """

        if self._optical_sigma is None:
            raise ValueError("optical_sigma is None")

        if self._voxel_size is None:
            raise ValueError("voxel_size not set")

        if self._pixel_size is None:
            raise ValueError("pixel_size not set")

        input = np.atleast_3d(np.array(input))
        if input.ndim > 3:
            raise TypeError("input can be 1d, 2d or 3d")

        shift = np.array(shift, int)
        if shift.size != 2 or shift.ndim != 1:
            raise TypeError("shift has to be (x,y)")

        if shift[0] > 0 and shift[1] > 0:
            input = input[shift[0]:, shift[1]:, ...]

        scale = self._voxel_size / self._pixel_size
        size = np.array(np.round(np.array(input.shape[0:2]) * scale), dtype=int)

        if np.amin(size) == 0:
            raise ValueError(
                f"voxel_size {self._voxel_size} and pixel_size {self._pixel_size} result in optical image size of {size}"
            )

        output = np.empty((size[0], size[1], input.shape[2]), dtype=input.dtype)

        if mp_pool and input.shape[2] > 1:
            chunk = [(input[:, :, i], self._optical_sigma, scale)
                     for i in range(input.shape[2])]

            results = mp_pool.starmap(optic.filter_resample, chunk)
            for i in range(input.shape[2]):
                output[:, :, i] = results[i]
        else:
            for i in range(input.shape[2]):
                output[:, :,
                       i] = optic.filter_resample(input[:, :, i],
                                                  self._optical_sigma, scale)

        return np.squeeze(output)

    def apply_untilt(self, input, theta, phi, mode='nearest'):
        """ applies optical untilt to image with affine transformation

        only necessary if untilt_view = False

        input: np.array(x,y(,rho))
        """

        if theta == 0:
            return input

        # calculate transformation matrix
        p = self._dim.copy()
        p_rot = 0.5 * np.array([p[0], p[1], 0])
        p_out = np.array([[p[0], p[1], 0], [p[0], 0, 0], [0, p[1], 0]])
        rot = tools.rotation.theta_phi(-theta, phi)

        p_in = np.array([np.dot(rot, p - p_rot) + p_rot for p in p_out])

        # TODO: refraction has to be implemented
        M = analysis.affine_transformation.calc_matrix(p_in[:, :2],
                                                       p_out[:, :2])

        output = analysis.affine_transformation.image(input, M, mode)
        if input.ndim == 3:
            output = np.atleast_3d(output)

        return output

    def apply_epa(self, input, mask=None):
        """ applies epa analysis to images

        input = np.array(x,y,rho)
        """

        transmittance, direction, retardation = analysis.epa.epa(input)
        if mask is not None:
            transmittance[np.invert(mask)] = float('nan')
            direction[np.invert(mask)] = float('nan')
            retardation[np.invert(mask)] = float('nan')

        return transmittance, direction, retardation

    def apply_rofl(self, input, mask=None, mp_pool=None, grad_mode=False):
        """ applies ROFL analysis to images

        input = np.array(tilt,x,y,rho)
        """

        self._print("Analyse tilts")

        if self._tilts is None:
            raise ValueError("tilts not set")

        if np.any(self._tilts[:, 1] != np.deg2rad([0, 0, 90, 180, 270])
                 ) or self._tilts[0, 0] != 0 or np.any(
                     self._tilts[1:, 0] != self._tilts[1, 0]):
            raise ValueError("tilts not suitable for ROFL")

        tilt_angle = self._tilts[1, 0]

        input = np.array(input, copy=False)

        if input.ndim != 4:
            raise TypeError("input: np.array([tilts,x,y,stack])")

        if input.shape[0] != 5:
            raise ValueError("input need 1 + 4 measurements")

        if input.shape[-1] <= 3:
            raise ValueError("input needs at least 3 equidistant rotations")

        if self._sensor_gain is None:
            raise ValueError("sensor_gain not defined")

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
                chunk = [(input[:, i, j, :], tilt_angle, self._sensor_gain, 0,
                          grad_mode) for i in range(input.shape[1])]
                results = mp_pool.starmap(analysis.rofl.rofl, chunk)

                for i, result in enumerate(results):
                    directionmap[i, j], inclmap[i, j], trelmap[i, j], dirdevmap[
                        i, j], incldevmap[i, j], treldevmap[i, j], funcmap[
                            i, j], itermap[i, j] = result
        else:
            for i in range(input.shape[1]):
                for j in range(input.shape[2]):
                    directionmap[i, j], inclmap[i, j], trelmap[i, j], dirdevmap[
                        i, j], incldevmap[i, j], treldevmap[i, j], funcmap[
                            i, j], itermap[i, j] = analysis.rofl.rofl(
                                input[:, i, j, :], tilt_angle,
                                self._sensor_gain, 0, grad_mode)

        return directionmap, inclmap, trelmap, (dirdevmap, incldevmap,
                                                treldevmap, funcmap, itermap)
