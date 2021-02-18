import fastpli.model.sandbox
import fastpli.model.solver
import fastpli.simulation
import fastpli.analysis
import fastpli.tools
import fastpli.io

import time
import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
import scipy.interpolate as si

pool = mp.Pool(2)
np.random.seed(42)

# two crossing main bundle
fb_0 = np.array([[-2000, -2000, 0], [-1000, -1000, 0], [0, 0, 0],
                 [1000, 1000, 0], [2000, 2000, 0]])
fb_1 = np.array([[-2000, 2000, 0], [-1000, 1000, 0], [0, 0, 0],
                 [1000, -1000, 0], [2000, -2000, 0]])

# left and right non crossing bundle
fb_left = np.array([[-2000, -2000, 0], [-1050, -1000, 0], [-300, 0, 0],
                    [-1050, 1000, 0], [-2000, 2000, 0]])
fb_right = fb_left.copy()
fb_right[:, 0] *= -1

# interpolate left and right bundle for smoother
N = 20
t = np.linspace(0, 1, fb_left.shape[0])
t_intp = np.linspace(0, 1, N)
fb_left = np.array([
    si.interp1d(t, fb_left[:, 0], 'quadratic')(t_intp),
    si.interp1d(t, fb_left[:, 1], 'quadratic')(t_intp),
    si.interp1d(t, fb_left[:, 2], 'quadratic')(t_intp)
]).T
fb_right = np.array([
    si.interp1d(t, fb_right[:, 0], 'quadratic')(t_intp),
    si.interp1d(t, fb_right[:, 1], 'quadratic')(t_intp),
    si.interp1d(t, fb_right[:, 2], 'quadratic')(t_intp)
]).T

fig, ax = plt.subplots(1, 1)
plt.plot(fb_0[:, 0], fb_0[:, 1])
plt.plot(fb_1[:, 0], fb_1[:, 1])
plt.plot(fb_left[:, 0], fb_left[:, 1])
plt.plot(fb_right[:, 0], fb_right[:, 1])
ax.axis('equal')
plt.show()

population = fastpli.model.sandbox.seeds.triangular_circle(800, 200)
fiber_radii = np.random.uniform(50.0, 100.0, population.shape[0])
ffb_0 = fastpli.model.sandbox.build.bundle(fb_0, population, fiber_radii)
ffb_1 = fastpli.model.sandbox.build.bundle(fb_1, population, fiber_radii)

population = fastpli.model.sandbox.seeds.triangular_circle(400, 200)
fiber_radii = np.random.uniform(50.0, 100.0, population.shape[0])
ffb_left = fastpli.model.sandbox.build.bundle(fb_left, population, fiber_radii)
ffb_right = fastpli.model.sandbox.build.bundle(fb_right, population,
                                               fiber_radii)

# Solving
solver = fastpli.model.solver.Solver()
solver.fiber_bundles = [ffb_0, ffb_1, ffb_left, ffb_right]
solver.obj_min_radius = 150
solver.obj_mean_length = 150
solver.omp_num_threads = 2

for i in range(1000):
    solved = solver.step()
    if i % (42 // 4) == 0:
        solver.draw_scene()
    if solved:
        break
solver.draw_scene()
time.sleep(4.2)
solver.close_scene()

# Simulation: generate discretised volume
simpli = fastpli.simulation.Simpli()
simpli.omp_num_threads = 2
simpli.voxel_size = 5  # in micro meter
simpli.set_voi([-2000, -2000, -50], [2000, 2000, 50])  # in micro meter
simpli.fiber_bundles = solver.fiber_bundles
simpli.fiber_bundles.layers = [[(1.0, -0.0001, 5, 'p')]] * len(
    simpli.fiber_bundles)
print('Memory:', str(round(simpli.memory_usage('MB'), 2)) + ' MB')

tissue, optical_axis, tissue_properties = simpli.generate_tissue()

# Simulation: measurment
simpli.filter_rotations = np.deg2rad([0, 30, 60, 90, 120, 150])
simpli.light_intensity = 26000  # a.u.
simpli.wavelength = 525  # in nm
simpli.pixel_size = 20  # in micro meter
simpli.optical_sigma = 0.71  # in pixel size
simpli.noise_model = lambda x: np.random.negative_binomial(x / (3 - 1), 1 / 3)
simpli.tilts = np.deg2rad([(0, 0), (5.5, 0), (5.5, 90), (5.5, 180), (5.5, 270)])

tilting_stack = []
for theta, phi in simpli.tilts:
    print(f'Tilt: theta={np.rad2deg(theta):.1f}, phi={np.rad2deg(phi):.1f}')
    images = simpli.run_simulation(tissue, optical_axis, tissue_properties,
                                   theta, phi)

    _, images = simpli.apply_optic(images)
    tilting_stack.append(images)

print('Generating modalities')
transmittance, direction, retardation = simpli.apply_epa(tilting_stack[0])

print('Tiltin analysis')
rofl_direction, rofl_incl, rofl_t_rel, _ = simpli.apply_rofl(tilting_stack,
                                                             mp_pool=pool)
fom = fastpli.analysis.images.fom_hsv_black(rofl_direction, rofl_incl)

# plot results
# Note: all data are xyz sorted.
fig, axs = plt.subplots(1, 3, sharey=True, figsize=(15, 5))
axs[0].imshow(transmittance.T, origin='lower', interpolation='nearest')
axs[0].set_title('transmittance')
axs[1].imshow(direction.T, origin='lower', interpolation='nearest')
axs[1].set_title('direction')
axs[2].imshow(retardation.T, origin='lower', interpolation='nearest')
axs[2].set_title('retardation')

fig, axs = plt.subplots(1, 4, sharey=True, figsize=(20, 5))
axs[0].imshow(rofl_direction.T, origin='lower', interpolation='nearest')
axs[0].set_title('direction')
axs[1].imshow(rofl_incl.T, origin='lower', interpolation='nearest')
axs[1].set_title('inclination')
axs[2].imshow(rofl_t_rel.T, origin='lower', interpolation='nearest')
axs[2].set_title('t_rel')
axs[3].imshow(np.swapaxes(fom, 0, 1), origin='lower', interpolation='nearest')
axs[3].set_title('fom')

fig, axs = plt.subplots(1, 1, figsize=(1, 1))
plt.imshow(np.swapaxes(fastpli.analysis.images.hsv_black_sphere(), 0, 1),
           origin='lower')
plt.title('fom - legend')
plt.axis('off')
plt.show()
