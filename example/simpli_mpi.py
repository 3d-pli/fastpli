import os
import time

import fastpli.simpli as simpli
import h5py
import numpy as np
from mpi4py import MPI

FILE_PATH = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # Read fiber data and prepair for PliGenerator
    # TODO: write json -> parameter function

    global_start = time.time()
    # construct templates
    pli_fiber_bundle = simpli.generation.FiberBundle()
    pli_fiber_layer_prop = [simpli.generation.LayerProperty(
        0.333, 0.004, 10, 1)]
    pli_fiber_layer_prop.append(simpli.generation.LayerProperty(
        0.666, -0.004, 5, 0))
    pli_fiber_layer_prop.append(simpli.generation.LayerProperty(
        1.0, 0.004, 1, 2))

    # read file
    with h5py.File(os.path.join(FILE_PATH, 'cube.h5'), 'r', driver='mpio', comm=MPI.COMM_WORLD) as h5f:
        fbs = h5f['fiber_bundles']
        for fb in fbs:
            fb = fbs[fb]
            for f in fb:
                f = fb[f]

                pli_fiber_data = simpli.generation.FiberData(
                    f['points'][:].flatten(), f['radii'][:])
                pli_fiber_bundle.push_fiber(pli_fiber_data)

        pli_fiber_bundle.set_fiber_bundle_properties(pli_fiber_layer_prop)

    # manipulation of fibers
    rot_mat = simpli.rotation.rot_zyz(
        np.deg2rad(20), np.deg2rad(-10), np.deg2rad(5))
    pli_fiber_bundle.rotate_fiber_bundle_around_point(
        list(rot_mat), [10, -5, 7.5])
    pli_fiber_bundle.translate_fiber_bundle([25, -15, 50])

    pli_fiber_bundles = [pli_fiber_bundle]

    # PliGenerator
    pli_gen = simpli.generation.Generator()
    dim = [300, 300, 300]
    pixel_size = 0.1

    global_dim = dim.copy()
    dim, pli_fiber_bundles, offset, overlap = simpli.mpi.get_local(
        dim, pli_fiber_bundles, np.deg2rad(5.5))
    offset += overlap

    pli_gen.set_volume(dim, pixel_size)
    pli_gen.set_fiber_bundle(pli_fiber_bundles)

    with h5py.File(os.path.join(FILE_PATH, 'output_mpi_' + str(comm.Get_size()) + '.h5'), 'w', driver='mpio', comm=MPI.COMM_WORLD)as h5f:
        h5f.atomic = True

        start = time.time()
        (label_field, vec_field,
         tissue_properties) = pli_gen.run_generation(0, 0)
        end = time.time()

        label_field_vis = label_field.asarray().astype(
            np.uint16).reshape(dim[0], dim[1], dim[2])
        label_field_vis[label_field_vis > 0] += 3
        dset = h5f.create_dataset(
            'tissue', (global_dim[0], global_dim[1], global_dim[2]), dtype=np.uint16)
        data = simpli.mpi.crop_data(label_field_vis, overlap)
        dset[offset[0]: offset[0] + data.shape[0],
             offset[1]: offset[1] + data.shape[1], :] = data

        dset = h5f.create_dataset(
            'vectorfield', (global_dim[0], global_dim[1], global_dim[2], 3), dtype=np.float32)
        data = simpli.mpi.crop_data(vec_field.asarray().reshape(
            dim[0], dim[1], dim[2], 3), overlap)
        dset[offset[0]: offset[0] + data.shape[0],
             offset[1]: offset[1] + data.shape[1], :] = data

        # PliSimulation
        pli_setup = simpli.simulation.PliSetup()
        pli_setup.filter_rotations = np.deg2rad([0, 30, 60, 90, 120, 150])
        pli_setup.light_intensity = 26000
        pli_setup.resolution = 1
        pli_setup.untilt_sensor = True
        pli_setup.wavelength = 525

        pli_sim = simpli.simulation.PliSimulator()
        pli_sim.set_pli_setup(pli_setup)
        pli_sim.set_tissue(label_field, vec_field, dim,
                           tissue_properties, pixel_size)

        if rank == 0:
            print("run_simulation: 0")
        image = pli_sim.run_simulation(0, 0, 1)
        image = image.reshape(dim[0], dim[1], len(pli_setup.filter_rotations))
        dset = h5f.create_dataset(
            'data/0', (global_dim[0], global_dim[1], image.shape[2]), dtype=np.float32)
        data = simpli.mpi.crop_data(image, overlap)

        dset[offset[0]:offset[0] + data.shape[0],
             offset[1]:offset[1] + data.shape[1]] = data

        if rank == 0:
            print("run_simulation: 1")
        image = pli_sim.run_simulation(np.deg2rad(5.5), np.deg2rad(0), 1)
        image = image.reshape(dim[0], dim[1], len(pli_setup.filter_rotations))
        dset = h5f.create_dataset(
            'data/1', (global_dim[0], global_dim[1], image.shape[2]), dtype=np.float32)
        data = simpli.mpi.crop_data(image, overlap)
        dset[offset[0]:offset[0] + data.shape[0],
             offset[1]:offset[1] + data.shape[1]] = data

        if rank == 0:
            print("run_simulation: 2")
        image = pli_sim.run_simulation(np.deg2rad(5.5), np.deg2rad(90), 1)
        image = image.reshape(dim[0], dim[1], len(pli_setup.filter_rotations))
        dset = h5f.create_dataset(
            'data/2', (global_dim[0], global_dim[1], image.shape[2]), dtype=np.float32)
        data = simpli.mpi.crop_data(image, overlap)
        dset[offset[0]:offset[0] + data.shape[0],
             offset[1]:offset[1] + data.shape[1]] = data

        if rank == 0:
            print("run_simulation: 3")
        image = pli_sim.run_simulation(np.deg2rad(5.5), np.deg2rad(180), 1)
        image = image.reshape(dim[0], dim[1], len(pli_setup.filter_rotations))
        dset = h5f.create_dataset(
            'data/3', (global_dim[0], global_dim[1], image.shape[2]), dtype=np.float32)
        data = simpli.mpi.crop_data(image, overlap)
        dset[offset[0]:offset[0] + data.shape[0],
             offset[1]:offset[1] + data.shape[1]] = data

        if rank == 0:
            print("run_simulation: 4")
        image = pli_sim.run_simulation(np.deg2rad(5.5), np.deg2rad(270), 1)
        image = image.reshape(dim[0], dim[1], len(pli_setup.filter_rotations))
        dset = h5f.create_dataset(
            'data/4', (global_dim[0], global_dim[1], image.shape[2]), dtype=np.float32)
        data = simpli.mpi.crop_data(image, overlap)
        dset[offset[0]:offset[0] + data.shape[0],
             offset[1]:offset[1] + data.shape[1]] = data

        end = time.time()
        # print("run_simulation:", end - start)

    global_end = time.time()

    if rank == 0:
        print("GlobalRuntime:", global_end - global_start)

    # print(max(image))
    # tracker = IndexTracker(ax, np.array(image).reshape(dim[0], dim[1], len(pli_setup.filter_rotations)))
    # fig.canvas.mpl_connect('scroll_event', tracker.onscroll)
    # plt.show()
