import fastpli.model.sandbox
import fastpli.model.solver
import fastpli.simulation
import fastpli.analysis
import fastpli.tools
import fastpli.io

import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
import imageio
import h5py
import os

pool = mp.Pool(2)
np.random.seed(42)

FILE_NAME = os.path.abspath(__file__)
FILE_PATH = os.path.dirname(FILE_NAME)
FILE_BASE = os.path.basename(FILE_NAME)
FILE_OUT = os.path.join(FILE_PATH, f'fastpli.example.{FILE_BASE}')


def data2image(data):
    return np.swapaxes(np.flip(data, 1), 0, 1)


fiber_bundle_trj_0 = [[-2000, 0, 0], [2000, 0, 0]]
fiber_bundle_trj_1 = [[0, -2000, 0], [0, 2000, 0]]

population = fastpli.model.sandbox.seeds.triangular_circle(750, 30)
population = fastpli.model.sandbox.seeds.crop_rectangle([-50, -10000],
                                                        [50, 10000],
                                                        population)

fiber_radii = np.random.uniform(10.0, 30.0, population.shape[0])
fiber_bundle_0 = fastpli.model.sandbox.build.bundle(fiber_bundle_trj_0,
                                                    population, fiber_radii)

population[:, [0, 1]] = population[:, [1, 0]]
fiber_radii = np.random.uniform(10.0, 30.0, population.shape[0])
fiber_bundle_1 = fastpli.model.sandbox.build.bundle(fiber_bundle_trj_1,
                                                    population, fiber_radii)

solver = fastpli.model.solver.Solver()
solver.fiber_bundles = [fiber_bundle_0, fiber_bundle_1]
solver.obj_min_radius = 8 * 20
solver.obj_mean_length = 4 * 20
solver.omp_num_threads = 2

# run solver
print('Begin solving process:')
solver.draw_scene()
for i in range(1000):
    solved = solver.step()

    # calculate current overlap ratio
    overlap = solver.overlap / solver.num_col_obj if solver.num_col_obj else 0
    if i % 10 == 0:
        print(f'step: {i}, {solver.num_obj}/{solver.num_col_obj} ' +
              f'{round(overlap * 100)}%')
        solver.draw_scene()
        # solver.save_ppm(f'solver_{i:03}.ppm')  # save a ppm image

    if i == 20:
        solver.draw_scene()
        solver.fiber_bundles = fastpli.objects.fiber_bundles.apply_fun_to_position(  # noqa: E501
            solver.fiber_bundles,
            lambda p: p + np.random.uniform(-10, 10, p.shape))
        solver.fiber_bundles = fastpli.objects.fiber_bundles.apply_fun_to_radii(  # noqa: E501
            solver.fiber_bundles,
            lambda r: r * np.random.lognormal(0, 0.1, r.shape))
        solver.draw_scene()

    if solved:
        print(f'solved: {i}, {solver.num_obj}/{solver.num_col_obj}')
        solver.draw_scene()
        break

fastpli.io.fiber_bundles.save(f'{FILE_OUT}.dat',
                              solver.fiber_bundles,
                              mode='w')
solver.close_scene()

# Simulation
print(f'creating file: {FILE_OUT}.h5')
with h5py.File(f'{FILE_OUT}.h5', 'w') as h5f:
    # save script
    h5f['version'] = fastpli.__version__
    with open(os.path.abspath(__file__), 'r') as f:
        h5f['parameter/script'] = f.read()
        h5f['parameter/pip_freeze'] = fastpli.tools.helper.pip_freeze()

    # Setup Simpli for Tissue Generation
    simpli = fastpli.simulation.Simpli()
    simpli.omp_num_threads = 2
    simpli.voxel_size = 2  # in micro meter
    simpli.set_voi([-2000, -2000, -30], [2000, 2000, 30])  # in micro meter
    simpli.fiber_bundles = fastpli.io.fiber_bundles.load(f'{FILE_OUT}.dat')

    # define layers (e.g. axon, myelin) inside fibers of each fiber_bundle
    simpli.fiber_bundles_properties = [[(1.0, -0.004, 10, 'p')]] * len(
        simpli.fiber_bundles)
    # (_0, _1, _2, _3)
    # _0: layer_scale times radius
    # _1: strength of birefringence
    # _2: absorption coefficient I = I*exp(-mu*x)
    # _3: model: 'p'-parallel, 'r'-radial or 'b'-background

    print('VOI:', simpli.get_voi())
    print('Memory:', str(round(simpli.memory_usage('MB'), 2)) + ' MB')

    # Generate Tissue
    print('Run Generation:')
    tissue, optical_axis, tissue_properties = simpli.generate_tissue()

    h5f['tissue/tissue'] = tissue.astype(np.uint16)
    h5f['tissue/optical_axis'] = optical_axis
    h5f['tissue/tissue_properties'] = tissue_properties

    with imageio.get_writer(f'{FILE_OUT}.tissue.gif', mode='I') as writer:
        m = np.amax(tissue)
        for i in range(tissue.shape[2]):
            writer.append_data(
                data2image(tissue[:, :, i] / m * 255).astype(np.uint8))

    # Simulate PLI Measurement
    simpli.filter_rotations = np.deg2rad(np.linspace(0, 180, 9,
                                                     endpoint=False))
    simpli.light_intensity = 26000  # a.u.
    simpli.interpolate = "Slerp"
    simpli.wavelength = 525  # in nm
    simpli.pixel_size = 20  # in micro meter
    simpli.optical_sigma = 0.71  # in pixel size
    simpli.noise_model = lambda x: np.random.negative_binomial(
        x / (3 - 1), 1 / 3)
    simpli.tilts = np.deg2rad([(0, 0), (5.5, 0), (5.5, 90), (5.5, 180),
                               (5.5, 270)])

    tilting_stack = [None] * 5
    print('Run Simulation:')
    for t, (theta, phi) in enumerate(simpli.tilts):
        print(round(np.rad2deg(theta), 1), round(np.rad2deg(phi), 1))
        images = simpli.run_simulation(tissue, optical_axis, tissue_properties,
                                       theta, phi)

        h5f['simulation/data/' + str(t)] = images

        # apply optic to simulation
        _, images = simpli.apply_optic(images)
        h5f['simulation/optic/' + str(t)] = images

        with imageio.get_writer(f'{FILE_OUT}.{t}.gif', mode='I') as writer:
            for i in range(images.shape[-1]):
                writer.append_data(
                    (data2image(images[:, :, i] - np.amin(images)) /
                     np.amax(images - np.amin(images)) * 255).astype(np.uint8))

        # calculate modalities
        epa = simpli.apply_epa(images)
        h5f[f'analysis/epa/{t}/transmittance'] = epa[0]
        h5f[f'analysis/epa/{t}/direction'] = np.rad2deg(epa[1])
        h5f[f'analysis/epa/{t}/retardation'] = epa[2]

        imageio.imwrite(f'{FILE_OUT}.{t}.transmittance.png',
                        data2image(epa[0]))
        imageio.imwrite(f'{FILE_OUT}.{t}.direction.png', data2image(epa[1]))
        imageio.imwrite(f'{FILE_OUT}.{t}.retardation.png', data2image(epa[2]))

        tilting_stack[t] = images

    # save mask for analysis
    mask = np.sum(tissue, 2) > 0
    mask = simpli.apply_optic_resample(1.0 * mask) > 0.1
    h5f['simulation/optic/mask'] = np.uint8(mask)
    mask = None  # keep analysing all pixels

    print('Run ROFL analysis:')
    rofl_direction, rofl_incl, rofl_t_rel, _ = simpli.apply_rofl(tilting_stack,
                                                                 mask=mask,
                                                                 mp_pool=pool)

    h5f['analysis/rofl/direction'] = np.rad2deg(rofl_direction)
    h5f['analysis/rofl/inclination'] = np.rad2deg(rofl_incl)
    h5f['analysis/rofl/trel'] = rofl_t_rel

    imageio.imwrite(f'{FILE_OUT}.r.direction.png', data2image(rofl_direction))
    imageio.imwrite(f'{FILE_OUT}.r.inclination.png', data2image(rofl_incl))
    imageio.imwrite(f'{FILE_OUT}.r.trel.png', data2image(rofl_t_rel))

    print(f'creating Fiber Orientation Map: {FILE_OUT}.png')

    imageio.imwrite(
        f'{FILE_OUT}.fom.png',
        data2image(
            fastpli.analysis.images.fom_hsv_black(rofl_direction, rofl_incl)))

print('simulation Done')
print('You can look at the data e.g. with Fiji and the hdf5 plugin')

# Plot
print('Plotting input and output orientations')
fig, axs = plt.subplots(1, 2, subplot_kw=dict(projection="polar"))
fbs_phi, fbs_theta = fastpli.analysis.orientation.fiber_bundles(
    solver.fiber_bundles)

pc = fastpli.analysis.orientation.histogram(fbs_phi,
                                            fbs_theta,
                                            axs[0],
                                            fun=lambda x: np.log(x + 1))
plt.colorbar(pc, ax=axs[0])

pc = fastpli.analysis.orientation.histogram(rofl_direction,
                                            np.pi / 2 - rofl_incl,
                                            axs[1],
                                            fun=lambda x: np.log(x + 1))
plt.colorbar(pc, ax=axs[1])

for ax in axs:
    ax.set_rmax(90)
    ax.set_rticks(range(0, 90, 10))
    ax.set_rlabel_position(22.5)
    ax.set_yticklabels([])
    ax.grid(True)
plt.show()
