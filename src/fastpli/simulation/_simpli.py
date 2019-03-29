# from . import __generation as generation
from .__generation import _Generator, _LayerProperty
from .__simulation import _Simulator, _Setup

from fastpli.objects import Fiber, Cell
from fastpli.tools import rotation

import numpy as np
import h5py
from mpi4py import MPI

# TODO: write json -> parameter function


class Simpli:

    def __init__(self):

        self._gen = _Generator()
        self._sim = _Simulator()
        self._fiber_bundles = None
        self._fiber_bundles_properties = None
        self._cells_populations = None
        self._cells_populations_properties = None
        self._dim = None
        self._dim_origin = np.array([0, 0, 0], dtype=int)
        self._pixel_size = None
        self._voi = None
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
        debug = bool(debug)
        self._debug = debug

    @property
    def dim(self):
        return self._dim

    @dim.setter
    def dim(self, dim):
        self._dim = np.array(dim)

    @property
    def dim_origin(self):
        return self._dim_origin

    @dim_origin.setter
    def dim_origin(self, dim_origin):
        self._dim_origin = np.array(dim_origin)

    @property
    def pixel_size(self):
        return (self._pixel_size)

    @pixel_size.setter
    def pixel_size(self, pixel_size):
        if not isinstance(pixel_size, (int, float)):
            raise TypeError("pixel_size : (int, float)")

        if pixel_size <= 0:
            raise ValueError("pixel_size !> 0")

        self._pixel_size = pixel_size

        if self._voi is not None:
            self.voi = self._voi
            self.print_debug("dim and dim_origin recalculated")

    @property
    def voi(self):
        if self._voi is None:
            if self._pixel_size is None:
                return None
                self.print_debug(
                    "pixel_size is not set, voi can't be calculated")

            if self._dim is None:
                return None
                self.print_debug("dim is not set, voi can't be calculated")

            self._voi = np.zeros((6,))
            self._voi[::2] = self._dim_origin + self._dim * self._pixel_size
            self._voi[1::2] = self._voi[::2] + self._dim * self._pixel_size

        return self._voi

    @voi.setter
    def voi(self, voi):
        voi = np.array(voi)
        if voi.size != 6 or voi.shape[0] != 6:
            raise TypeError("voi: wrong shape, has to be (6,)")

        if voi[0] > voi[1] or voi[2] > voi[3] or voi[4] > voi[5]:
            raise ValueError(
                "voi not corrected sorted: (x_min, x_max, y_min, y_max, z_min, z_max)"
            )

        self._voi = voi

        if self._pixel_size is None:
            self.print_debug(
                "pixel_size is not set yet, dim and dim_origin will be calculated when setted"
            )
            return

        tmp = np.array(self._voi / self._pixel_size)
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
                if isinstance(f, Fiber):
                    continue
                elif isinstance(f, (list, tuple)):
                    f = np.array(f)

                if isinstance(f, np.ndarray):
                    if len(f.shape) is not 2 or f.shape[1] is not 4:
                        raise TypeError("fiber elements has to be of dim nx4")
                    fbs[fb_i][f_i] = Fiber(f[:, 0:-1], f[:, -1])
                else:
                    raise TypeError(
                        "fiber hast to be a objects.Fiber, 4d-list or 4d-array")

        self._fiber_bundles = fbs

    @property
    def fiber_bundles_properties(self):
        return self._fiber_bundles_properties

    @fiber_bundles_properties.setter
    def fiber_bundles_properties(self, bundle_layer_properties):
        if not isinstance(bundle_layer_properties, (list, tuple)):
            raise TypeError("properties must be a list(list(tuples))")

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
                        Fiber(f['points'][:].flatten(), f['radii'][:]))

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
        self._gen.set_volume(self._dim, self.dim_origin, self.pixel_size)
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
        setup = _Setup()
        setup.light_intensity = self._light_intensity
        setup.pixel_size = self._pixel_size
        setup.wavelength = self._wavelength
        setup.untilt_sensor = self._untilt_sensor
        setup.filter_rotations = self._filter_rotations
        self._sim.set_pli_setup(setup)

    def RunSimulation(self,
                      label_field,
                      vec_field,
                      tissue_properties,
                      theta,
                      phi,
                      do_untilt=True):

        self.InitSimulation()
        self._sim.set_tissue_properties(tissue_properties)
        image = self._sim.run_simulation(self._dim, label_field, vec_field,
                                         theta, phi, do_untilt)
        return image

    def DimData(self):
        dim_local = self._gen.dim_local()
        dim_offset = self._gen.dim_offset()
        return dim_local, dim_offset

    def SaveAsH5(self, h5f, data, data_name):

        dim_local, dim_offset = self.DimData()
        if data_name is 'tissue':
            dim = self._dim
            dset = h5f.create_dataset(data_name, dim, dtype=np.uint16)

            for i in range(data.shape[0]):
                dset[i + dim_offset[0], dim_offset[1]:dim_offset[1] +
                     dim_local[1], dim_offset[2]:dim_offset[2] +
                     dim_local[2]] = data[i, :, :]

        elif data_name is 'vectorfield':
            dim = [self._dim[0], self._dim[1], self._dim[2], 3]
            dset = h5f.create_dataset(data_name, dim, dtype=np.float32)
            for i in range(data.shape[0]):
                dset[i + dim_offset[0], dim_offset[1]:dim_offset[1] +
                     dim_local[1], dim_offset[2]:dim_offset[2] +
                     dim_local[2]] = data[i, :, :]
        elif 'data/' in data_name:
            dim = [self._dim[0], self._dim[1], self._filter_rotations.size]
            dset = h5f.create_dataset(data_name, dim, dtype=np.float32)

            if tuple(dim) == data.shape:
                dset[:] = data[:]
            else:
                mask = (np.count_nonzero(data, axis=2) != 0)

                for i in range(data.shape[0]):

                    first = 0
                    for idx, elm in enumerate(mask[i, :]):
                        if elm:
                            first = idx
                            break

                    last = -1
                    for idx, elm in reversed(list(enumerate(mask[i, :]))):
                        if elm:
                            last = idx + 1
                            break

                    if first <= last:
                        dset[i + dim_offset[0], first + dim_offset[1]:last +
                             dim_offset[1], :] = data[i, first:last, :]

        else:
            raise TypeError("no compatible SaveAsH5: " + data_name)
