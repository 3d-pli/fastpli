# -*- coding: utf-8 -*-
"""
Simpli CPU Class
"""

from .__generation import _Generator
from .__simulation import _Simulator

from . import _mpi
from .. import analysis

import numpy as np
import warnings

from .__simpli import __Simpli


class Simpli(__Simpli):
    """
    Simpli Class for simulating 3D-PLI images
    """

    def __init__(self, mpi_comm=None):

        super().__init__()

        # LIBRARIES
        self.__gen = _Generator()
        self.__sim = _Simulator()
        self.mpi = None
        if mpi_comm:
            from mpi4py import MPI
            self.__gen.set_mpi_comm(MPI._addressof(mpi_comm))
            self.__sim.set_mpi_comm(MPI._addressof(mpi_comm))
            self.mpi = _mpi._MPI(mpi_comm)

        # freeze class
        self.__freeze()

    @property
    def interpolate(self):
        """ (de)activate interpolation of vectors in simulation: bool """
        return self._interpolate

    @interpolate.setter
    def interpolate(self, interpolate):
        if interpolate != 'NN' and interpolate != 'Lerp' and \
           interpolate != 'Slerp':
            raise ValueError('Only \'NN\', \'Lerp\' or \'Slerp\' are supported')
        self._interpolate = interpolate

    def generate_tissue(self, only_tissue=False):
        """ generating discret tissue for simulation """

        self._print('Generate Tissue')
        self._print(f'Memory needed: ~{np.ceil(self.memory_usage()):.0f} MB')

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
            only_tissue)

        if not np.any(tissue):
            warnings.warn(
                'All labels are 0. Usually this means, that the VOI contains no\
                     fiber_bundles or that the voxel_size is to large.',
                UserWarning)

        if self.mpi:
            self.mpi._set_gen_dim(self.__gen.dim_local(), self._dim,
                                  self.__gen.dim_offset())

        return tissue, optical_axis, tissue_properties

    def _init_pli_setup(self):
        """ check all pli setup input requirements """

        self._check_simulation_input()
        self.__sim.set_pli_setup(self._step_size, self._light_intensity,
                                 self._voxel_size, self._wavelength,
                                 self._tissue_refrection, self._interpolate,
                                 self._untilt_sensor_view, self._flip_z_beam,
                                 self._filter_rotations)

    def _check_generation_input(self):
        if not self._fiber_bundles and not self._cells_populations:
            raise ValueError('fiber_bundles and cells_populations are not set')
        self._check_property_length()

    def run_simulation(self, tissue, optical_axis, tissue_properties, theta,
                       phi):
        """
        running simulation. Input from tissue_generation

        tissue and optical_axis will be passed by reference
        theta, phi: tilting angle in radiant
        """

        self._print('Run simulation')
        self._check_volume_input()
        self._init_pli_setup()

        tissue_properties = np.array(tissue_properties)
        if tissue_properties.ndim != 2:
            raise ValueError('tissue_properties.ndim !=2')

        if tissue_properties.shape[1] != 2:
            raise ValueError('tissue_properties.shape[1] != 2')

        tissue_ = np.array(tissue, dtype=np.int32, copy=False)
        optical_axis_ = np.array(optical_axis, dtype=np.float32, copy=False)

        if tissue_ is not tissue:
            warnings.warn('tissue is copied', UserWarning)
        if optical_axis_ is not optical_axis:
            warnings.warn('optical_axis is copied', UserWarning)

        if tissue.ndim != 3:
            raise ValueError('tissue.ndim != 3')
        if optical_axis.ndim != 4:
            raise ValueError('optical_axis.ndim != 4')
        if optical_axis.shape[3] != 3:
            raise ValueError('optical_axis.shape[3] != 3')
        if not np.all(tissue.shape == optical_axis.shape[:-1]):
            raise ValueError(
                'not np.equal(tissue.shape, optical_axis.shape[:-1])')
        if np.any(tissue.shape != self._dim):
            raise ValueError('np.any(tissue.shape != self._dim)')

        images = self.__sim.run_simulation(self._dim, tissue_, optical_axis_,
                                           tissue_properties, theta, phi)
        if np.min(images.flatten()) < 0:
            raise ValueError('intensity < 0 detected')

        return images

    def run_tissue_pipeline(self,
                            h5f=None,
                            script=None,
                            save=['tissue', 'optical_axis']):
        """ Automatic pipeline for tissue generation with save options """

        self._print('Run tissue pipeline')
        self._check_volume_input()
        self._check_generation_input()

        if script:
            self.save_parameter_h5(h5f, script)

        if 'all' in save or 'simulation' in save:
            save = save + ['tissue', 'optical_axis']

        # run tissue generation
        tissue, optical_axis, tissue_properties = self.generate_tissue()

        if h5f and 'tissue' in save:
            self._print('Save tissue')
            dset = h5f.create_dataset('tissue/tissue',
                                      tissue.shape,
                                      dtype=np.uint16,
                                      compression='gzip',
                                      compression_opts=1)
            dset[:] = tissue

        if h5f and 'optical_axis' in save:
            self._print('Save optical_axis')
            dset = h5f.create_dataset('tissue/optical_axis',
                                      optical_axis.shape,
                                      dtype=np.float32,
                                      compression='gzip',
                                      compression_opts=1)
            dset[:] = optical_axis

        if h5f and 'tissue' in save:
            h5f['tissue/properties'] = tissue_properties

        return tissue, optical_axis, tissue_properties

    def crop_tilt_pixel(self):
        """ crop affected boundary pixel from tilted images """
        self._print('calc crop-tilt-halo pixel')
        if not self.untilt_sensor_view:
            raise ValueError('currently only for untilt_sensor_view=True')
        if self._voxel_size is None:
            raise ValueError('voxel_size is not set')
        if self._pixel_size is None:
            raise ValueError('pixel_size is not set')
        if self._tilts is None:
            raise ValueError('tilts are not set')

        delta_pixel = int(
            np.ceil(
                np.tan(np.max(np.abs(self._tilts[:, 0]))) * self._dim[-1] / 2 *
                self._voxel_size / self._pixel_size))
        return delta_pixel

    def crop_tilt_voxel(self):
        """ get number of affected boundary voxel from tilted images """
        self._print('calc crop-tilt-halo voxel')
        delta_voxel = int(
            np.round(self.crop_tilt_pixel() * self._pixel_size /
                     self._voxel_size))
        return delta_voxel

    def add_crop_tilt_halo(self):
        """ add number of necessary boundary voxel from tilted images """
        self._print('add crop-tilt-halo')
        self.dim_origin[:2] -= self.crop_tilt_voxel() * self.voxel_size
        self.dim[:2] += 2 * self.crop_tilt_voxel()

    def rm_crop_tilt_halo(self, input):
        """ remove number of added boundary voxel from tilted images """
        self._print('rm crop-tilt-halo')
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
        """
        Automatic pipeline for simulation and analysis with save options
        """

        self._print('Run simulation pipeline')
        if 'all' in save or 'simulation' in save:
            save = save + ['data', 'resample', 'optic', 'epa', 'mask', 'rofl']

        if self._tilts is None:
            raise ValueError('tilts is not set')
        if self._optical_sigma is None:
            raise ValueError('optical_sigma is not set')

        flag_rofl = True
        if np.any(self._tilts[:, 1] != np.deg2rad([0, 0, 90, 180, 270])
                 ) or self._tilts[0, 0] != 0 or np.any(
                     self._tilts[1:, 0] != self._tilts[1, 0]):
            warnings.warn('Tilts not suitable for ROFL. Skipping analysis')
            flag_rofl = False

        tilting_stack = [None] * len(self._tilts)

        self._print('Simulate tilts:')
        for t, tilt in enumerate(self._tilts):
            theta, phi = tilt[0], tilt[1]
            self._print(f'Tilt {t}: theta: {np.rad2deg(theta):.1f} deg, phi: ' +
                        f'{np.rad2deg(phi):.1f} deg')
            images = self.run_simulation(tissue, optical_axis,
                                         tissue_properties, theta, phi)

            if crop_tilt:
                images = self.rm_crop_tilt_halo(images)

            if h5f and 'data' in save:
                self._print('Save data')
                h5f['simulation/data/' + str(t)] = images
                h5f['simulation/data/' + str(t)].attrs['theta'] = theta
                h5f['simulation/data/' + str(t)].attrs['phi'] = phi

            # apply optic to simulation
            img_res, img_noise = self.apply_optic(images, mp_pool=mp_pool)

            if h5f and 'resample' in save:
                self._print('Save resample')
                h5f['simulation/resample/' + str(t)] = img_res
                h5f['simulation/resample/' + str(t)].attrs['theta'] = theta
                h5f['simulation/resample/' + str(t)].attrs['phi'] = phi

            if h5f and 'optic' in save:
                self._print('Save optic')
                h5f['simulation/optic/' + str(t)] = img_noise
                h5f['simulation/optic/' + str(t)].attrs['theta'] = theta
                h5f['simulation/optic/' + str(t)].attrs['phi'] = phi

            if h5f and 'epa' in save:
                # calculate modalities
                epa = self.apply_epa(img_noise)

                self._print('Save epa')
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

            tilting_stack[t] = img_noise

        # pseudo mask
        if h5f and 'mask' in save:
            self._print('Calc mask')
            mask = np.sum(tissue, 2) > 0
            mask = self.apply_optic_resample(1.0 * mask, mp_pool=mp_pool) > 0.1
            self._print('Save mask')
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
            self._print('Save rofl')
            h5f['analysis/rofl/direction'] = rofl_direction
            h5f['analysis/rofl/inclination'] = rofl_incl
            h5f['analysis/rofl/t_rel'] = rofl_t_rel

        if h5f and flag_rofl and 'rofl_conf' in save:
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

    def run_pipeline(self,
                     h5f=None,
                     script=None,
                     save=[
                         'tissue', 'optical_axis', 'data', 'optic', 'epa',
                         'mask', 'rofl'
                     ],
                     crop_tilt=False,
                     mp_pool=None):
        """
        Automatic tissue generation and simulation pipeline with save options
        """

        self._print('Run pipeline')
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
            raise ValueError('tilts is not set')
        if self._optical_sigma is None:
            raise ValueError('optical_sigma is not set')

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
            raise TypeError('num_threads != int')

        if num_threads <= 0:
            raise TypeError('num_threads <= 0')

        num_threads_gen = self.__gen.set_omp_num_threads(num_threads)
        num_threads_sim = self.__sim.set_omp_num_threads(num_threads)

        if num_threads_gen != num_threads_sim:
            raise AssertionError('num_threads_gen != num_threads_sim')

        if num_threads_gen != num_threads:
            warnings.warn('reduced num_threads: ' + str(num_threads_gen),
                          UserWarning)

        self._omp_num_threads = num_threads_gen

    def memory_usage(self, unit='MB', item='all'):
        """ print expected memory ussage """
        if unit == 'MB':
            div = 1024**2
        elif unit == 'GB':
            div = 1024**3
        else:
            raise ValueError('allowed is only \'MB\', \'GB\'')

        if self._dim is None:
            raise TypeError('dimension not set yet')

        if item == 'tissue':
            # tissue + distance_array
            return np.prod(self._dim) * (32 + 32) / 8 / div
        elif item == 'all':
            return np.prod(self._dim) * (32 + 32 + 3 * 32) / 8 / div
        else:
            raise ValueError('allowed is only \'tissue\' or \'all\'')
