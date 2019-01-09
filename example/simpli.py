import time
import h5py
import fastpli.simpli as simpli
import numpy as np
import sys
import os

# FILE_PATH = '/home/homeGlobal/fmatuschke/code/fastPLI/example'
FILE_PATH = os.path.dirname(os.path.abspath(__file__))
# sys.path.insert(0, FILE_PATH+'/../build')


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
with h5py.File(os.path.join(FILE_PATH, 'cube.h5'), 'r') as h5f:
    fbs = h5f['fiber_bundles']
    for fb in fbs:
        fb = fbs[fb]
        for f in fb:
            f = fb[f]

            pli_fiber_data = simpli.generation.FiberData(
                f['points'][:].flatten(), f['radii'][:])
            pli_fiber_bundle.push_fiber(pli_fiber_data)

    pli_fiber_bundle.set_fiber_bundle_properties(pli_fiber_layer_prop)

with h5py.File(os.path.join(FILE_PATH, 'output.h5'), 'w') as h5f:

    # manipulation of fibers
    rot_mat = simpli.rotation.rot_zyz(
        np.deg2rad(20), np.deg2rad(-10), np.deg2rad(5))
    pli_fiber_bundle.rotate_fiber_bundle_around_point(
        list(rot_mat), [10, -5, 7.5])
    pli_fiber_bundle.translate_fiber_bundle([25, -15, 50])

    # PliGenerator
    pli_gen = simpli.generation.Generator()

    dim = [100, 100, 100]
    pixel_size = 1
    pli_gen.set_volume(dim, pixel_size)
    pli_gen.set_fiber_bundle([pli_fiber_bundle])

    start = time.time()
    (label_field, vec_field,
     tissue_properties) = pli_gen.run_generation(0, 0)
    end = time.time()

    tmp = label_field.asarray().astype(
        np.uint16).reshape(dim[0], dim[1], dim[2])
    tmp[tmp > 0] += 3
    h5f['tissue'] = tmp

    dset = h5f.create_dataset(
        'vectorfield', (dim[0], dim[1], dim[2], 3), dtype=np.float32)
    simpli.data_container.save_2_h5(dset, vec_field, dim + [3])

    print("TissueGeneration:", end - start)

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

    start = time.time()
    print("run_simulation: 0")
    image = pli_sim.run_simulation(0, 0, 1)
    h5f['data/0'] = image.reshape(dim[0],
                                  dim[1], len(pli_setup.filter_rotations))

    print("run_simulation: 1")
    image = pli_sim.run_simulation(np.deg2rad(5.5), np.deg2rad(0), 1)
    h5f['data/1'] = image.reshape(dim[0],
                                  dim[1], len(pli_setup.filter_rotations))

    print("run_simulation: 2")
    image = pli_sim.run_simulation(np.deg2rad(5.5), np.deg2rad(90), 1)
    h5f['data/2'] = image.reshape(dim[0],
                                  dim[1], len(pli_setup.filter_rotations))

    print("run_simulation: 3")
    image = pli_sim.run_simulation(np.deg2rad(5.5), np.deg2rad(180), 1)
    h5f['data/3'] = image.reshape(dim[0],
                                  dim[1], len(pli_setup.filter_rotations))

    print("run_simulation: 4")
    image = pli_sim.run_simulation(np.deg2rad(5.5), np.deg2rad(270), 1)
    h5f['data/4'] = image.reshape(dim[0],
                                  dim[1], len(pli_setup.filter_rotations))

    end = time.time()
    print("run_simulation:", end - start)

global_end = time.time()
print("GlobalRuntime:", global_end - global_start)

# print(max(image))
# tracker = IndexTracker(ax, np.array(image).reshape(dim[0], dim[1], len(pli_setup.filter_rotations)))
# fig.canvas.mpl_connect('scroll_event', tracker.onscroll)
# plt.show()
