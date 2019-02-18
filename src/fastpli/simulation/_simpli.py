#!/usr/bin/python3

from fastpli.simulation import generation, simulation
from fastpli.objects import Fiber
from fastpli.tools import rotation
import numpy as np
import h5py

# Read fiber data and prepair for PliGenerator
# TODO: write json -> parameter function


class Simpli:

    def __init__(self):
        self.__fiber_bundles = [[]]
        self.__fiber_bundles_properties = [[]]

        self.__gen = generation.Generator()
        self.__sim = simulation.Simulator()

        # self.layer_properties = [[]]
        self.dim = [0, 0, 0]
        self.dim_origin = [0, 0, 0]

        self.pixel_size = 0
        self.setup = simulation.Setup()

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
            self.dim, self.dim_origin, self.pixel_size)
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

        image = self.__sim.run_simulation(self.dim,
                                          label_field, vec_field, theta, phi, do_untilt)
        return image

    def DimData(self):
        dim_local = self.__gen.dim_local()
        dim_offset = self.__gen.dim_offset()
        return dim_local, dim_offset
