import fastpli.model.sandbox
import fastpli.io
import fastpli.tools
import fastpli.objects
import fastpli.model

import numpy as np
from tqdm import tqdm, trange

### 3 crossing fiber bundles ###

RADIUS_OUT = 10000
RADIUS_IN = RADIUS_OUT * 0.6
FIBER_SPACING = 25
FIBER_RADIUS = 10
FIBER_STEPS = 100

p0 = np.array((0, 0, -20))
p1 = np.array((0, 0, 80))
p_shift = np.array(
    ((RADIUS_OUT + RADIUS_IN) * 0.5 / np.cos(np.deg2rad(30)), 0, 0))

fiber_bundles = []
### first ###
dp = np.dot(fastpli.tools.rotation.z(np.deg2rad(-120)), p_shift)
data = fastpli.model.sandbox.shape.cylinder(p0 + dp, p1 + dp, RADIUS_IN,
                                            RADIUS_OUT, np.deg2rad(30),
                                            np.deg2rad(30 + 60), 'c',
                                            FIBER_SPACING, FIBER_STEPS)

fiber_bundles.append(
    fastpli.objects.add_radius_to_fiber_bundle(data, FIBER_RADIUS).copy())

### second ###
dp = np.dot(fastpli.tools.rotation.z(np.deg2rad(120)), p_shift)
data = fastpli.model.sandbox.shape.cylinder(p0 + dp, p1 + dp, RADIUS_IN,
                                            RADIUS_OUT, np.deg2rad(-30 - 60),
                                            np.deg2rad(-30), 'c', FIBER_SPACING,
                                            FIBER_STEPS)
fiber_bundles.append(
    fastpli.objects.add_radius_to_fiber_bundle(data, FIBER_RADIUS).copy())

### third ###
dp = p_shift
data = fastpli.model.sandbox.shape.cylinder(p0 + dp, p1 + dp, RADIUS_IN,
                                            RADIUS_OUT, np.deg2rad(150),
                                            np.deg2rad(150 + 60), 'c',
                                            FIBER_SPACING, FIBER_STEPS)
fiber_bundles.append(
    fastpli.objects.add_radius_to_fiber_bundle(data, FIBER_RADIUS).copy())

for fb in fiber_bundles:
    for f in fb:
        f[:, :-1] += np.random.uniform(-FIBER_RADIUS * 0.5, FIBER_RADIUS * 0.5,
                                       (f.shape[0], 3))

### SS ###
dp = p_shift
data = fastpli.model.sandbox.shape.cylinder(p0 + dp, p1 + dp, RADIUS_IN * 0.6,
                                            RADIUS_IN, np.deg2rad(150),
                                            np.deg2rad(150 + 60), 'p',
                                            FIBER_SPACING, 2)
fiber_bundles.append(
    fastpli.objects.add_radius_to_fiber_bundle(data, FIBER_RADIUS).copy())

### radial ###
dp = np.dot(fastpli.tools.rotation.z(np.deg2rad(-120)), p_shift)
data = fastpli.model.sandbox.shape.cylinder(p0 + dp, p1 + dp, RADIUS_IN * 0.6,
                                            RADIUS_IN * 1.2, np.deg2rad(30),
                                            np.deg2rad(30 + 60), 'r',
                                            FIBER_SPACING, 20)
fiber_bundles.append(
    fastpli.objects.add_radius_to_fiber_bundle(data, FIBER_RADIUS).copy())

for fb in fiber_bundles:
    for f in fb:
        f[:, :-1] += np.random.uniform(-FIBER_RADIUS * 0.5, FIBER_RADIUS * 0.5,
                                       (f.shape[0], 3))

fastpli.io.fiber.save('test.dat', fiber_bundles)

solver = fastpli.model.Solver()
solver.omp_num_threads = 8
solver.fiber_bundles = fiber_bundles
solver.obj_mean_length = FIBER_RADIUS * 2
solver.obj_min_radius = FIBER_RADIUS * 4
for i in trange(10000):
    solved = solver.step()
    if solved:
        break

    if (i % 100) == 0:
        print("step:", i, solver.num_obj, solver.num_col_obj)
        fastpli.io.fiber.save('test_' + str(i) + '.dat', solver.fiber_bundles)
        # solver.draw_scene()

fastpli.io.fiber.save('test_.dat', solver.fiber_bundles)
