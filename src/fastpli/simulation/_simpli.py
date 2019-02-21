#!/usr/bin/python3

from fastpli.simulation import generation, simulation
from fastpli.objects import Fiber
from fastpli.tools import rotation
import numpy as np
import h5py
from mpi4py import MPI

# Read fiber data and prepair for PliGenerator
# TODO: write json -> parameter function


class Simpli:

    def __init__(self):
        self.__fiber_bundles = [[]]
        self.__fiber_bundles_properties = [[]]

        self.__gen = generation.Generator()
        self.__sim = simulation.Simulator()

        # self.layer_properties = [[]]
        self._dim = None
        self._dim_origin = np.array([0, 0, 0], dtype=int)

        self.pixel_size = 0
        self.setup = simulation.Setup()

    @property
    def dim(self):
        return self._dim
    
    @dim.setter
    def dim(self,dim):
        self._dim = np.array(dim)

    @property
    def dim_origin(self):
        return self._dim_origin
    
    @dim_origin.setter
    def dim_origin(self,dim_origin):
        self._dim_origin = np.array(dim_origin)

    def ReadFiberFile(self, filename):
        self.__fiber_bundles = []
        with h5py.File(filename, 'r') as h5f:

            fbs = h5f['fiber_bundles']
            for fb in fbs:
                self.__fiber_bundles.append([])
                fb = fbs[fb]
                for f in fb:
                    f = fb[f]
                    self.__fiber_bundles[-1].append(
                        Fiber(f['points'][:].flatten(), f['radii'][:]))

    def TranslateVolume(self, offset):
        if not isinstance(offset, list):
            raise TypeError("offset must be a list")

        for fb in self.__fiber_bundles:
            for f in fb:
                f.translate(offset)

    def RotateVolume(self, phi, theta, psi):
        rot_mat = rotation.zyz(phi, theta, psi)
        for fb in self.__fiber_bundles:
            for f in fb:
                f.rotate(list(rot_mat))

    def RotateVolumeAroundPoint(self, phi, theta, psi, offset):
        offset = np.array(offset)
        if offset.shape != (3,):
            raise TypeError("offset must a point")
        rot_mat = rotation.zyz(phi, theta, psi)
        for fb in self.__fiber_bundles:
            for f in fb:
                f.rotate_around_point((rot_mat), offset)

    def SetFiberProperties(self, bundle_layer_properties):
        if not isinstance(bundle_layer_properties, list):
            raise TypeError("properties must be a list(list(tuples))")

        if len(self.__fiber_bundles) != len(bundle_layer_properties):
            raise TypeError(
                "properties must have the same size as fiber_bundles")

        self.__fiber_bundles_properties = []
        for prop in bundle_layer_properties:
            if not isinstance(prop, list):
                raise TypeError(
                    "properties must be a list(list(tuples))")

            self.__fiber_bundles_properties.append([])
            for ly in prop:
                if len(ly) != 4:
                    raise TypeError(
                        "properties must have len 4 (float, float, float, char)")
                self.__fiber_bundles_properties[
                    -1].append(generation.LayerProperty(ly[0], ly[1], ly[2], ly[3]))

    def __CheckFiberBundleAndPropertiesLength(self):
        if len(self.__fiber_bundles) != len(self.__fiber_bundles_properties):
            raise TypeError(
                "properties must have the same size as fiber_bundles")

    def GenerateTissue(self, only_label=False):
        self.__gen.set_volume(
            self._dim, self.dim_origin, self.pixel_size)
        self.__CheckFiberBundleAndPropertiesLength()
        self.__gen.set_fiber_bundles(
            self.__fiber_bundles, self.__fiber_bundles_properties)
        label_field, vec_field, tissue_properties = self.__gen.run_generation(
            only_label, 0)

        return label_field, vec_field, tissue_properties

    def InitSimulation(self, label_field, vec_field, tissue_properties):
        self.__sim.set_pli_setup(self.setup)
        self.__sim.set_tissue(
            label_field, vec_field, tissue_properties, self.pixel_size)

    def RunSimulation(self, label_field, vec_field,
                      tissue_properties, theta, phi, do_untilt=True):

        self.setup.pixel_size = self.pixel_size
        self.__sim.set_pli_setup(self.setup)
        self.__sim.set_tissue_properties(tissue_properties)

        image = self.__sim.run_simulation(self._dim,
                                          label_field, vec_field, theta, phi, do_untilt)
        return image

    def DimData(self):
        dim_local = self.__gen.dim_local()
        dim_offset = self.__gen.dim_offset()
        return dim_local, dim_offset

    def SaveAsH5(self, h5f, data, data_name):

        dim_local, dim_offset = self.DimData()
        if data_name is 'tissue':
            dim = self._dim
            dset = h5f.create_dataset(data_name, dim, dtype=np.uint16)

            for i in range(data.shape[0]):
                dset[i+dim_offset[0], dim_offset[1]:dim_offset[1]
                 + dim_local[1], dim_offset[2]:dim_offset[2] + dim_local[2]] = data[i,:,:]

        elif data_name is 'vectorfield':
            dim = [self._dim[0], self._dim[1], self._dim[2], 3]
            dset = h5f.create_dataset(data_name, dim, dtype=np.float32)
            for i in range(data.shape[0]):
                dset[i+dim_offset[0], dim_offset[1]:dim_offset[1]
                 + dim_local[1], dim_offset[2]:dim_offset[2] + dim_local[2]] = data[i,:,:]
        elif 'data/' in data_name:
            dim = [self._dim[0], self._dim[1], len(self.setup.filter_rotations)]
            dset = h5f.create_dataset(data_name, dim, dtype=np.float32)

            if tuple(dim) == data.shape:
                dset[:] = data[:]
            else:
                mask = (np.count_nonzero(data, axis=2) != 0)

                for i in range(data.shape[0]):

                    first = 0
                    for idx, elm in enumerate(mask[i,:]):
                        if elm:
                            first = idx
                            break
                    

                    last = -1
                    for idx, elm in reversed(list(enumerate(mask[i,:]))):
                        if elm:
                            last = idx+1
                            break
                    
                    if first <= last:
                        dset[i+dim_offset[0], first+dim_offset[1]:last+dim_offset[1],:] = data[i, first:last,:]
                

        else:
            raise TypeError("no compatible SaveAsH5: " + data_name)
