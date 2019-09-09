import os
import numpy as np
import h5py


def load(file_name, dataset_name='/'):
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
    elif ext == '.h5':
        with h5py.File(file_name, 'r') as h5f:
            if dataset_name[-1] is not '/':
                dataset_name = dataset_name + '/'
            fbs = h5f[dataset_name]
            for i in range(len(fbs)):
                if i != 0:
                    fiber_bundles.append([])
                fb = fbs[str(i)]
                for j in range(len(fb)):
                    fiber_bundles[-1].append(fb[str(j) +
                                                '/data'][:].astype(float))
    else:
        raise TypeError(ext + ' is not implemented yet')

    return fiber_bundles


def save(file_name, fiber_bundles, dataset_name='/'):
    _, ext = os.path.splitext(file_name)

    if ext == '.dat' or ext == '.txt':
        with open(file_name, 'w') as file:
            for fb, fiber_bundle in enumerate(fiber_bundles):
                for fiber in fiber_bundle:

                    np.array(fiber, dtype=np.float64, copy=False)

                    if fiber.shape[1] != 4 or len(fiber.shape) != 2:
                        raise TypeError('Wrong shape:', fiber.shape)

                    for line in fiber:
                        file.write(
                            str(line[0]) + " " + str(line[1]) + " " +
                            str(line[2]) + " " + str(line[3]) + "\n")
                    file.write("\n")
                if fb != len(fiber_bundles) - 1:
                    file.write("\n")
    elif ext == '.h5':
        with h5py.File(file_name, 'w') as h5f:
            if dataset_name[-1] is not '/':
                dataset_name = dataset_name + '/'

            for fb_i, fb in enumerate(fiber_bundles):
                for i, f in enumerate(fb):
                    h5f[dataset_name + str(fb_i) + '/' + str(i) +
                        '/data'] = f[:, :]
    else:
        raise TypeError(ext + ' is not implemented yet')
