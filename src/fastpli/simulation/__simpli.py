# -*- coding: utf-8 -*-
"""
Simpli Class
"""

from . import optic
from .. import analysis
from .. import objects
from .. import tools
from .. import __version__
from .. import __compiler__
from .. import __libraries__

import numpy as np
import warnings
import sys


class __Simpli:
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
        self._interpolate = "Slerp"
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
        """ origin of volume in micro: (3)-array """
        return self._dim_origin

    @dim_origin.setter
    def dim_origin(self, dim_origin):
        self._dim_origin = np.array(dim_origin, dtype=float)

    @property
    def voxel_size(self):
        """ size of voxel in micro: float """
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
        """ pixel size of resulting optical image in micro: float """
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

        if sensor_gain < 0:
            raise ValueError("sensor_gain is < 0")

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

        if optical_sigma < 0:
            raise ValueError("optical_sigma is < 0")

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

                if ly[1] < 0 and ly[-1] == 'r':
                    warnings.warn("birefringence negative and radial")
                if ly[1] > 0 and ly[-1] == 'p':
                    warnings.warn("birefringence positive and parallel")
                if ly[1] != 0 and ly[-1] == 'b':
                    warnings.warn(
                        "birefringence != 0 for background. Will be set to 0")

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

    def save_parameter_h5(self, h5f, script=None):
        """ Saves class members without fiber_bundles in hdf5 file. """
        self._print("Save fastpli parameter")
        h5f.attrs['fastpli/simpli'] = str(self.get_dict())
        h5f.attrs['fastpli/version'] = __version__
        h5f.attrs['fastpli/compiler'] = __compiler__
        h5f.attrs['fastpli/libraries'] = __libraries__
        h5f.attrs['fastpli/pip_freeze'] = tools.helper.pip_freeze()
        h5f.attrs['fastpli/system'] = sys.version
        if script:
            h5f.attrs['script'] = script.encode().decode('ascii', 'ignore')

    def apply_optic(self, input, shift=(0, 0), mp_pool=None):
        """ applies optical resample, convolution and noise to image
        input: np.array(x,y(,rho))
        """

        self._print("Apply optic")
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

        self._print("Apply optic resample")
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
                f"voxel_size {self._voxel_size}, pixel_size {self._pixel_size} \
                    and input shape {input.shape[0:2]} result in optical image \
                    size of {size}")

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

        self._print("Apply untilt")
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

        self._print("Apply epa")
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

        self._print("Apply rofl")
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
                    if not mask[i, j]:
                        continue
                    directionmap[i, j], inclmap[i, j], trelmap[i, j], dirdevmap[
                        i, j], incldevmap[i, j], treldevmap[i, j], funcmap[
                            i, j], itermap[i, j] = result
        else:
            for i in range(input.shape[1]):
                for j in range(input.shape[2]):
                    if not mask[i, j]:
                        continue
                    directionmap[i, j], inclmap[i, j], trelmap[i, j], dirdevmap[
                        i, j], incldevmap[i, j], treldevmap[i, j], funcmap[
                            i, j], itermap[i, j] = analysis.rofl.rofl(
                                input[:, i, j, :], tilt_angle,
                                self._sensor_gain, 0, grad_mode)

        return directionmap, inclmap, trelmap, (dirdevmap, incldevmap,
                                                treldevmap, funcmap, itermap)
