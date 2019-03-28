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
                    fiber.append(numbers)
                if not line.strip():
                    if fiber:
                        fiber_bundles[-1].append(np.array(fiber))
                    fiber = []
    else:
        raise TypeError(ext + ' is not implemented yet')

    return fiber_bundles
