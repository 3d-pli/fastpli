# -*- coding: utf-8 -*-
"""
Simpli Class
"""

import warnings
import sys

import numpy as np

from . import optic
from .. import analysis
from .. import objects
from .. import tools
from .. import __version__
from .. import __compiler__
from .. import __libraries__


class __Simpli:
    """
    Simpli Class for simulating 3D-PLI images
    """

    __is_frozen = False

    def __setattr__(self, key, value):
        if self.__is_frozen and not hasattr(self, key):
            raise TypeError(f'{self.__class__.__name__} is a frozen class')
        object.__setattr__(self, key, value)

    def __freeze(self):
        self.__is_frozen = True

    def __init__(self):

        # DIMENSIONS
        self._dim = None
        self._dim_origin = np.array([0, 0, 0], dtype=float)
        self._voxel_size = None

        # GENERATION
        self._fiber_bundles = None

        # SIMULATION
        self._filter_rotations = None
        self._flip_z_beam = False
        self._interpolate = 'Slerp'
        self._light_intensity = None
        self._pixel_size = None
        self._optical_sigma = None
        self._step_size = 1.0
        self._tilts = None
        self._tissue_refrection = 1
        self._background_absorption = 0
        self._untilt_sensor_view = True
        self._wavelength = None

        # ANALYSIS
        self._noise_model = None

        # OTHER
        self._omp_num_threads = 1
        self._verbose = False

    def _print(self, msg):
        if self._verbose:
            print(msg)

    def get_dict(self):
        """ Get all member variables which are properties as dictionary """
        members = dict()
        for key, value in self.__dict__.items():
            if key in ('_fiber_bundles'):
                continue
            if key.startswith('_') and not key.startswith(
                    '__') and not key.startswith('_Simpli'):

                if isinstance(value, np.ndarray):
                    members[key[1:]] = value.tolist()
                else:
                    members[key[1:]] = value

        return members

    def set_dict(self, data):
        """ Set dictionary of variables to class members """
        for key, value in data.items():
            if key.startswith('_'):
                raise ValueError('member variable cant be set directly')

            if value is not None:
                setattr(self, key, value)
            else:
                warnings.warn('None value in dict detected')

    @property
    def dim(self):
        """ dimension of volume in voxel: (3)-np.ndarray """
        return self._dim

    @dim.setter
    def dim(self, dim):
        dim = np.array(dim)
        if dim.dtype != int:
            raise TypeError('dim is not np.ndarray(int)')
        if dim.size != 3:
            raise TypeError('dim.size != 3')
        if dim.ndim != 1:
            raise TypeError('dim.ndim != 1')

        self._dim = dim

    @property
    def dim_origin(self):
        """ origin of volume in micro meter: (3)-np.ndarray """
        return self._dim_origin

    @dim_origin.setter
    def dim_origin(self, dim_origin):
        self._dim_origin = np.array(dim_origin, dtype=float)

    @property
    def voxel_size(self):
        """ size length of voxel in micro meter """
        return self._voxel_size

    @voxel_size.setter
    def voxel_size(self, voxel_size):
        if voxel_size <= 0:
            raise ValueError('voxel_size <= 0')

        flag = False
        if self._voxel_size is not None and self._dim is not None and \
           self._dim_origin is not None:
            v_min, v_max = self.get_voi()
            flag = True

        self._voxel_size = float(voxel_size)

        if flag and self._dim is not None and self._dim_origin is not None:
            self.set_voi(v_min, v_max)
            self._print('voxel_size: recalculated dimensions')

    @property
    def pixel_size(self):
        """ side length of pixel of resulting optical image in micro meter """
        return self._pixel_size

    @pixel_size.setter
    def pixel_size(self, pixel_size):
        if pixel_size <= 0:
            raise ValueError('pixel_size <= 0')
        self._pixel_size = float(pixel_size)

    def get_voi(self):
        """
        get volume of interest

        Returns
        -------
        res: np.ndarray, np.ndarray
            v_min: [x_min, y_min, z_min]
            v_max: [x_max, y_max, z_max]
        """
        if self._voxel_size is None:
            raise ValueError('voxel_size is not set, voi can\'t be calculated')

        if self._dim is None:
            raise ValueError('dim is not set, voi can\'t be calculated')

        if self._dim_origin is None:
            raise ValueError('dim_origin is not set, voi can\'t be calculated')

        voi = np.zeros((6,))
        voi[::2] = self._dim_origin
        voi[1::2] = voi[::2] + self._dim * self._voxel_size

        v_min = np.array(voi[0::2])
        v_max = np.array(voi[1::2])
        return v_min, v_max

    def set_voi(self, v_min, v_max):
        """
        set volume of interest

        Parameters
        ----------
        v_min: (float,float,float)-array-like
            [x_min, y_min, z_min]
        v_max: (float,float,float)-array-like
            [x_max, y_max, z_max]
        """

        v_min = np.array(v_min, dtype=float)
        v_max = np.array(v_max, dtype=float)

        if v_min.ndim != 1 or v_max.ndim != 1:
            raise TypeError('v_min,v_max : ndim != 1')
        if v_min.size != 3 or v_max.size != 3:
            raise TypeError('v_min,v_max : size != 1')
        if np.any(v_min >= v_max):
            raise ValueError('v_min >= v_max')

        if self._voxel_size is None:
            raise TypeError('voxel_size is not set yet')

        self._dim = np.array(np.round((v_max - v_min) / self._voxel_size),
                             dtype=int)
        self._dim_origin = v_min

    @property
    def filter_rotations(self):
        """ list of filter rotation position in radiant: [float] """
        return self._filter_rotations

    @filter_rotations.setter
    def filter_rotations(self, filter_rotations):
        filter_rotations = np.array(filter_rotations, dtype=float)

        if filter_rotations.size == 0 or filter_rotations.ndim != 1:
            raise TypeError('filter_rotations : nx1')

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
        """
        list of spherical tilting angles [[theta, phi], ...] in radiant:

        Parameters
        ----------
        tilts: [[float, float], ...]
        """
        return self._tilts

    @tilts.setter
    def tilts(self, tilts):
        if not isinstance(tilts, (list, tuple, np.ndarray)):
            raise TypeError('tilts is not a list or array')

        tilts = np.array(tilts, ndmin=2)
        if tilts.ndim != 2:
            raise TypeError('tilts.shape != nx2')
        if tilts.shape[-1] != 2:
            raise TypeError('tilts.shape != nx2')

        self._tilts = tilts

    @property
    def noise_model(self):
        """
        noise model to apply on resampled image

        Parameters
        ----------
        fun: function
        """
        return self._noise_model

    @noise_model.setter
    def noise_model(self, noise_model):

        if not callable(noise_model):
            raise TypeError('noise model is not callable')

        self._noise_model = noise_model

    @property
    def optical_sigma(self):
        """
        optical_sigma of convolution to applie before resampled image

        Parameters
        ----------
        optical_sigma: float
            in pixel size
        """
        return self._optical_sigma

    @optical_sigma.setter
    def optical_sigma(self, optical_sigma):

        if not isinstance(optical_sigma, (int, float)):
            raise TypeError('optical_sigma is not a number')

        if optical_sigma < 0:
            raise ValueError('optical_sigma is < 0')

        self._optical_sigma = optical_sigma

    @property
    def wavelength(self):
        """
        wavelength of light

        Parameters
        ----------
        wavelength: float
            in nano meter
        """
        return self._wavelength

    @wavelength.setter
    def wavelength(self, wavelength):
        self._wavelength = float(wavelength)

    @property
    def background_absorption(self):
        """
        absorption coefficient of the background

        Parameters
        ----------
        background_absorption: float
        """
        return self._background_absorption

    @background_absorption.setter
    def background_absorption(self, background_absorption):
        self._background_absorption = float(background_absorption)

    @property
    def tissue_refrection(self):
        """
        tissue refrection to include in tilting light path calculation

        Parameters
        ----------
        tissue_refrection: float
        """
        return self._tissue_refrection

    @tissue_refrection.setter
    def tissue_refrection(self, tissue_refrection):
        self._tissue_refrection = float(tissue_refrection)

    @property
    def step_size(self):
        """
        step_size of light inside tissue in voxeln

        Parameters
        ----------
        step_size: float
            of light inside tissue in voxeln
        """
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
        """
        automatic untilting of images

        Parameters
        ----------
        untilt_sensor_view: bool
            True or False
        """
        return self._untilt_sensor_view

    @untilt_sensor_view.setter
    def untilt_sensor_view(self, untilt_sensor_view):
        self._untilt_sensor_view = bool(untilt_sensor_view)

    @property
    def flip_z_beam(self):
        """
        flip light z-direction

        Parameters
        ----------
        flip_z_beam: bool
            True or False
        """
        return self._flip_z_beam

    @flip_z_beam.setter
    def flip_z_beam(self, flip_z_beam):
        self._flip_z_beam = bool(flip_z_beam)

    @property
    def fiber_bundles(self):
        """ get fiber_bundles: FiberBundles """
        return self._fiber_bundles

    @fiber_bundles.setter
    def fiber_bundles(self, fbs):
        """ set fiber_bundles: [[(,4)-array]] """
        if fbs is None:
            self._fiber_bundles = None
        else:
            fbs = objects.FiberBundles(fbs)
            self._fiber_bundles = fbs

    def _check_volume_input(self):
        if self._dim is None:
            raise ValueError('dim not set')

        if self._dim_origin is None:
            raise ValueError('dim_origin not set')

        if self._voxel_size is None:
            raise ValueError('voxel_size not set')

        self.get_voi()

    def _check_generation_data(self):
        if not self._fiber_bundles:
            raise ValueError('fiber_bundles is not set')

    def _check_simulation_input(self):
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

    def save_parameter_h5(self, h5f, script=None):
        """ Saves class members without fiber_bundles in hdf5 file. """
        self._print('Save fastpli parameter')
        h5f.attrs['fastpli/simpli'] = str(self.get_dict())
        h5f.attrs['fastpli/version'] = __version__
        h5f.attrs['fastpli/compiler'] = __compiler__
        h5f.attrs['fastpli/libraries'] = __libraries__
        h5f.attrs['fastpli/pip_freeze'] = tools.helper.pip_freeze()
        h5f.attrs['fastpli/system'] = sys.version
        if script:
            h5f.attrs['script'] = script.encode().decode('ascii', 'ignore')

    def apply_optic(self, data, mp_pool=None):
        """
        applies optical filters and resampling to image

        Parameters
        ----------
        data: np.ndarray
            images(x,y,rho)
        mp_pool: optional, multiprocessing.Pool

        Return
        ------
        res: np.ndarray
            images in reduced size and applied noise model
        """

        self._print('Apply optic')
        if self._noise_model is None:
            raise ValueError('noise_model not set')

        resampled = self.apply_optic_resample(data, mp_pool=mp_pool)

        if np.amin(resampled) < 0:
            raise AssertionError('intensity < 0 detected')

        noisy = optic.add_noise(resampled, self._noise_model)

        return resampled, noisy

    def apply_optic_resample(self, data, shift=(0, 0), mp_pool=None):
        """
        applies optical resampling to image

        Parameters
        ----------
        data: np.ndarray
            images(x,y,rho)
        mp_pool: optional, multiprocessing.Pool

        Return
        ------
        res: np.ndarray
            images in reduced size
        """

        self._print('Apply optic resample')
        if self._optical_sigma is None:
            raise ValueError('optical_sigma is None')

        if self._voxel_size is None:
            raise ValueError('voxel_size not set')

        if self._pixel_size is None:
            raise ValueError('pixel_size not set')

        data = np.atleast_3d(np.array(data))
        if data.ndim > 3:
            raise TypeError('data can be 1d, 2d or 3d')

        shift = np.array(shift, int)
        if shift.size != 2 or shift.ndim != 1:
            raise TypeError('shift has to be (x,y)')

        if shift[0] > 0 and shift[1] > 0:
            data = data[shift[0]:, shift[1]:, ...]

        scale = self._voxel_size / self._pixel_size
        size = np.array(np.round(np.array(data.shape[0:2]) * scale), dtype=int)

        if np.amin(size) == 0:
            raise ValueError(f'voxel_size {self._voxel_size}, ' +
                             f'pixel_size {self._pixel_size} ' +
                             f'and data shape {data.shape[0:2]} ' +
                             f'result in optical image size of {size}')

        output = np.empty((size[0], size[1], data.shape[2]), dtype=data.dtype)

        if mp_pool and data.shape[2] > 1:
            chunk = [(data[:, :, i], self._optical_sigma, scale)
                     for i in range(data.shape[2])]

            results = mp_pool.starmap(optic.filter_resample, chunk)
            for i in range(data.shape[2]):
                output[:, :, i] = results[i]
        else:
            for i in range(data.shape[2]):
                output[:, :,
                       i] = optic.filter_resample(data[:, :, i],
                                                  self._optical_sigma, scale)

        return np.squeeze(output)

    def apply_untilt(self, data, theta, phi, mode='nearest'):
        """
        applies optical untilt to image with affine transformation

        Parameters
        ----------
        data: np.ndarray
            images(x,y,rho)
        theta: float
            polar angle
        phi: float
            azimuthal angle
        mode: str
            method of scipy.interpolate.griddata

        Return
        ------
        res: np.ndarray
            untilted images
        """

        if theta == 0:
            return data

        self._print('Apply untilt')
        # calculate transformation matrix
        p = self._dim.copy()
        p_rot = 0.5 * np.array([p[0], p[1], 0])
        p_out = np.array([[p[0], p[1], 0], [p[0], 0, 0], [0, p[1], 0]])
        rot = tools.rotation.zymz(-theta, phi)

        p_in = np.array([np.dot(rot, p - p_rot) + p_rot for p in p_out])

        # TODO: refraction has to be implemented
        M = analysis.affine_transformation.calc_matrix(p_in[:, :2],
                                                       p_out[:, :2])

        output = analysis.affine_transformation.apply(data, M, mode)
        if data.ndim == 3:
            output = np.atleast_3d(output)

        return output

    def apply_epa(self, data, mask=None):
        """
        calculates transmittance, direction and retardation

        Parameters
        ----------
        data: np.ndarray
            images(x,y,rho)
        mask: np.ndarray(bool)
            only True elements will be calculated

        Return
        ------
        res: np.ndarray, np.ndarray, np.ndarray
            transmittance, direction, retardation
        """

        self._print('Apply epa')
        transmittance, direction, retardation = analysis.epa.epa(data)
        if mask is not None:
            transmittance[np.invert(mask)] = float('nan')
            direction[np.invert(mask)] = float('nan')
            retardation[np.invert(mask)] = float('nan')

        return transmittance, direction, retardation

    def apply_rofl(self, data, mask=None, mp_pool=None, grad_mode=False):
        """
        calculates tilt analysis for direction, inclination and t_rel and
        additional fitting parameters

        Parameters
        ----------
        data: np.ndarray
            images(x,y,rho)
        mask: optional, np.ndarray(bool)
            only True elements will be calculated
        mp_pool: optional, multiprocessing.Pool
        grad_mode: optional, bool
            activate gradient method inside rofl algorithm

        Return
        ------
        res: np.ndarray, np.ndarray, np.ndarray,
             (np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray)
            directionmap, inclmap, trelmap,
            (dirdevmap, incldevmap, treldevmap, funcmap, itermap)
        """

        self._print('Apply rofl')
        if self._tilts is None:
            raise ValueError('tilts not set')

        if np.any(self._tilts[:, 1] != np.deg2rad([0, 0, 90, 180, 270])
                 ) or self._tilts[0, 0] != 0 or np.any(
                     self._tilts[1:, 0] != self._tilts[1, 0]):
            raise ValueError('tilts not suitable for ROFL')

        tilt_angle = self._tilts[1, 0]

        data = np.array(data, copy=False)
        dtype = np.float32 if data.itemsize <= 4 else np.float64
        data = data.astype(dtype)

        if data.ndim != 4:
            raise TypeError('data: np.array([tilts,x,y,stack])')

        if data.shape[0] != 5:
            raise ValueError('data need 1 + 4 measurements')

        if data.shape[-1] <= 3:
            raise ValueError('data needs at least 3 equidistant rotations')

        if mask is None:
            mask = np.ones((data.shape[1], data.shape[2]), bool)

        dtype = np.float32 if data.itemsize <= 4 else np.float64
        directionmap = np.empty_like(mask, dtype=dtype)
        inclmap = np.empty_like(mask, dtype=dtype)
        trelmap = np.empty_like(mask, dtype=dtype)
        dirdevmap = np.empty_like(mask, dtype=dtype)
        incldevmap = np.empty_like(mask, dtype=dtype)
        treldevmap = np.empty_like(mask, dtype=dtype)
        funcmap = np.empty_like(mask, dtype=dtype)
        itermap = np.empty_like(mask, dtype=dtype)

        if mp_pool:
            for j in range(data.shape[2]):
                chunk = [(data[:, i, j, :], tilt_angle, 3, 0, grad_mode)
                         for i in range(data.shape[1])]
                results = mp_pool.starmap(analysis.rofl.rofl, chunk)

                for i, result in enumerate(results):
                    if not mask[i, j]:
                        continue
                    directionmap[i, j], inclmap[i, j], trelmap[i, j], dirdevmap[
                        i, j], incldevmap[i, j], treldevmap[i, j], funcmap[
                            i, j], itermap[i, j] = result
        else:
            for i in range(data.shape[1]):
                for j in range(data.shape[2]):
                    if not mask[i, j]:
                        continue
                    directionmap[i, j], inclmap[i, j], trelmap[i, j], dirdevmap[
                        i, j], incldevmap[i, j], treldevmap[i, j], funcmap[
                            i, j], itermap[i, j] = analysis.rofl.rofl(
                                data[:, i, j, :], tilt_angle, 3, 0, grad_mode)

        return directionmap, inclmap, trelmap, (dirdevmap, incldevmap,
                                                treldevmap, funcmap, itermap)
