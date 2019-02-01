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

    def run_simulation(self, label_field, vec_field,
                       tissue_properties, theta, phi, do_untilt=True):

        self.__sim.set_pli_setup(self.setup)
        self.__sim.set_tissue(
            self.dim, tissue_properties, self.pixel_size)

        image = self.__sim.run_simulation(
            label_field, vec_field, theta, phi, do_untilt)
        return image.reshape(self.dim[0], self.dim[1], len(self.setup.filter_rotations))
    
    def DimData(self):
        dim_local = self.__gen.dim_local()
        dim_offset = self.__gen.dim_offset()
        return dim_local, dim_offset


if __name__ == "__main__":

    with h5py.File('output.h5', 'w') as h5f:
        simpli = Simpli()

        # PliGeneration ###
        simpli.pixel_size = 1
        simpli.dim = [100, 100, 100]
        simpli.ReadFiberFile('example/simpli/cube.h5')
        simpli.layer_properties = [[0.333, 0.004, 10, 1], [
            0.666, -0.004, 5, 0], [1.0, 0.004, 1, 2]]

        # manipulation of fibers
        simpli.RotateVolumeAroundPoint(np.deg2rad(
            20), np.deg2rad(-10), np.deg2rad(5), [10, -5, 7.5])
        simpli.TranslateVolume([25, -15, 50])

        label_field, vec_field, tissue_properties = simpli.GenerateTissue()

        label_field_vis = np.transpose(label_field.asarray().astype(
            np.uint16).reshape(simpli.dim[0], simpli.dim[1], simpli.dim[2]), (0, 1, 2))
        label_field_vis[label_field_vis > 0] += 3
        h5f['tissue'] = label_field_vis
        h5f['vectorfield'] = np.transpose(vec_field.asarray().reshape(
            simpli.dim[0], simpli.dim[1], simpli.dim[2], 3), (0, 1, 2, 3))

        # PliSimulation ###
        simpli.setup.filter_rotations = np.deg2rad([0, 30, 60, 90, 120, 150])
        simpli.setup.light_intensity = 26000
        simpli.setup.resolution = 1
        simpli.setup.untilt_sensor = True
        simpli.setup.wavelength = 525

        simpli.InitSimulation(label_field, vec_field, tissue_properties)

        print("run_simulation: 0")
        image = simpli.run_simulation(0, 0)
        h5f['data/0'] = np.transpose(image, (0, 1, 2))

        print("run_simulation: 1")
        image = simpli.run_simulation(np.deg2rad(5.5), np.deg2rad(0))
        h5f['data/1'] = np.transpose(image, (0, 1, 2))

        print("run_simulation: 2")
        image = simpli.run_simulation(np.deg2rad(5.5), np.deg2rad(90))
        h5f['data/2'] = np.transpose(image, (0, 1, 2))

        print("run_simulation: 3")
        image = simpli.run_simulation(np.deg2rad(5.5), np.deg2rad(180))
        h5f['data/3'] = np.transpose(image, (0, 1, 2))

        print("run_simulation: 4")
        image = simpli.run_simulation(np.deg2rad(5.5), np.deg2rad(270))
        h5f['data/4'] = np.transpose(image, (0, 1, 2))
