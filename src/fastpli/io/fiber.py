import os
import numpy as np
from .. import objects


def load(file_name):
    _, ext = os.path.splitext(file_name)

    fiber_bundles = [[]]
    if ext == '.dat' or ext == '.txt':
        with open(file_name, 'r') as f:
            fiber = []
            flag_fiber_bundle_end = False
            for line in f:
                if line.strip():
                    if flag_fiber_bundle_end:
                        fiber_bundles.append([])
                        flag_fiber_bundle_end = False

                    numbers = list(map(float, line.split()))
                    fiber.append(numbers[0:4])
                if not line.strip():
                    if fiber:
                        fiber_bundles[-1].append(np.array(fiber))
                        fiber = []
                    else:  # new bundle with double empty line
                        flag_fiber_bundle_end = True
            if fiber:
                fiber_bundles[-1].append(np.array(fiber))
    else:
        raise TypeError(ext + ' is not implemented yet')

    return fiber_bundles


def save(file_name, fiber_bundles):
    _, ext = os.path.splitext(file_name)

    if ext == '.dat' or ext == '.txt':
        with open(file_name, 'w') as file:
            for fb, fiber_bundle in enumerate(fiber_bundles):
                for fiber in fiber_bundle:
                    if isinstance(fiber, objects.Fiber):
                        pos, radii = fiber.data
                    else:
                        if isinstance(fiber, list):
                            fiber = np.array(fiber)
                        if not isinstance(fiber, np.ndarray):
                            raise TypeError('Wrong input datatype')

                        if fiber.shape[1] != 4 or len(fiber.shape) != 2:
                            raise TypeError('Wrong shape:', fiber.shape)

                        pos = fiber[:, 0:3]
                        radii = fiber[:, -1]

                    for i in range(len(pos)):
                        file.write(
                            str(pos[i, 0]) + " " + str(pos[i, 1]) + " " +
                            str(pos[i, 2]) + " " + str(radii[i]) + "\n")
                    file.write("\n")
                if fb != len(fiber_bundles) - 1:
                    file.write("\n")

    else:
        raise TypeError(ext + ' is not implemented yet')
