#!/usr/bin/python3

from fastpli.simpli import generation, simulation, rotation
import numpy as np
import h5py

# Read fiber data and prepair for PliGenerator
# TODO: write json -> parameter function


class Simpli:

    def __init__(self):
        self.__fiber_bundles = []
        self.__gen = generation.Generator()
        self.__sim = simulation.PliSimulator()

        self.layer_properties = [[]]
        self.dim = [0, 0, 0]
        self.pixel_size = 0
        self.setup = simulation.PliSetup()

    def ReadFiberFile(self, filename):
        with h5py.File(filename, 'r') as h5f:
            fiber_bundle = generation.FiberBundle()

            fbs = h5f['fiber_bundles']
            for fb in fbs:
                fb = fbs[fb]
                for f in fb:
                    f = fb[f]
                    data = generation.FiberData(
                        f['points'][:].flatten(), f['radii'][:])
                    fiber_bundle.push_fiber(data)

                self.__fiber_bundles.append(fiber_bundle)

    def TranslateVolume(self, offset):
        if not isinstance(offset, list):
            raise TypeError("offset must be a list")

        for fb in self.__fiber_bundles:
            fb.translate_fiber_bundle(offset)

    def RotateVolume(self, phi, theta, psi):
        rot_mat = rotation.rot_zyz(phi, theta, psi)
        for fb in self.__fiber_bundles:
            fb.rotate_fiber_bundle(list(rot_mat))

    def RotateVolumeAroundPoint(self, phi, theta, psi, offset):
        if not isinstance(offset, list):
            raise TypeError("offset must be a list")
        rot_mat = rotation.rot_zyz(phi, theta, psi)
        for fb in self.__fiber_bundles:
            fb.rotate_fiber_bundle_around_point(list(rot_mat), offset)

    def __SetFiberProperties(self):
        if not isinstance(self.layer_properties, list):
            raise TypeError("properties must be a list(list)")
        if not isinstance(self.layer_properties[0], list):
            raise TypeError("properties must be a list(list)")

        if isinstance(self.layer_properties[0][0], list):
            if len(self.layer_properties) is not len(self.__gen.fiber_bundles):
                raise TypeError(
                    "propertie length and fiber_bundle length is not the same")
            for fb, pp in zip(self.__gen.fiber_bundles, self.layer_properties):
                lp = []
                for p in pp:
                    lp.append(generation.LayerProperty(p))
                fb.set_fiber_bundle_properties(lp)
        else:
            for fb in self.__fiber_bundles:
                lp = []
                for p in self.layer_properties:
                    lp.append(generation.LayerProperty(p))
                fb.set_fiber_bundle_properties(lp)

    def GenerateTissue(self):
        self.__gen.set_volume(self.dim, self.pixel_size)
        self.__SetFiberProperties()
        self.__gen.set_fiber_bundle(self.__fiber_bundles)
        label_field, vec_field, tissue_properties = self.__gen.run_generation(
            0, 0)

        return label_field, vec_field, tissue_properties

    def InitSimulation(self, label_field, vec_field, tissue_properties):
        print("InitSimualtion:", self.pixel_size, self.dim)
        self.__sim.set_pli_setup(self.setup)
        self.__sim.set_tissue(label_field, vec_field, self.dim,
                              tissue_properties, self.pixel_size)

    def run_simulation(self, theta, phi, do_untilt=True):

        image = self.__sim.run_simulation(theta, phi, do_untilt)
        return image.reshape(self.dim[0], self.dim[1],
                             len(self.setup.filter_rotations))


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
