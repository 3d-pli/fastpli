import os
import numpy as np


def read(file_name):
    _, ext = os.path.splitext(file_name)

    fiber_bundles = [[]]
    if ext == '.dat' or ext == '.txt':
        with open(file_name, 'r') as f:
            fiber = []
            for line in f:
                if line.strip():
                    numbers = list(map(float, line.split()))
                    fiber.append(numbers[0:4])
                if not line.strip():
                    if fiber:
                        fiber_bundles[-1].append(np.array(fiber))
                        flag_fiber_end = True
                    else:  # new bundle with double empty line
                        fiber_bundles.append([])
                    fiber = []
    else:
        raise TypeError(ext + ' is not implemented yet')

    return fiber_bundles


def save(file_name, fiber_bundles):
    _, ext = os.path.splitext(file_name)

    if ext == '.dat' or ext == '.txt':
        with open(file_name, 'w') as file:
            for fiber_bundle in fiber_bundles:
                for fiber in fiber_bundle:
                    pos, radii = fiber.data
                    for i in range(len(pos)):
                        file.write(
                            str(pos[i, 0]) + " " + str(pos[i, 1]) + " " +
                            str(pos[i, 2]) + " " + str(radii[i]) + "\n")
                    file.write("\n")
                file.write("\n")

    else:
        raise TypeError(ext + ' is not implemented yet')

    return fiber_bundles
