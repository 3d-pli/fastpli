import os
import numpy as np
import h5py


def load(file_name, group_name='fiber_bundles/'):
    """
    Load fiberbundles configurations from a text file oder hdf5 file

    Parameters
    ----------
    file_name : string
        path to file
    group_name : string, optional
        path inside hdf5-file for fiber bundles

    Returns
    -------
    res : list(list(fiber)), fibers are (n,4)-arrays with (x,y,z,radii) for each fiber point
    """
    _, ext = os.path.splitext(file_name)

    fiber_bundles = []
    if ext in ['.dat', '.txt']:
        with open(file_name, 'r') as f:
            fiber = []
            fiber_bundles.append([])
            for line in f:
                if line.strip():
                    numbers = list(map(float, line.split()))
                    fiber.append(numbers[0:4])
                if not line.strip():
                    if fiber:
                        fiber_bundles[-1].append(np.array(fiber, float))
                        fiber = []
                    else:  # new bundle
                        fiber_bundles.append([])
            if fiber:
                fiber_bundles[-1].append(np.array(fiber))
    elif ext == '.h5':
        with h5py.File(file_name, 'r') as h5f:
            fb_list = list(map(int, list(h5f[group_name])))
            fb_list.sort()
            for fb in fb_list:
                fiber_bundles.append([])
                f_list = h5f[group_name][str(fb)]
                f_list = list(map(int, list(h5f[group_name][str(fb)])))
                f_list.sort()
                for f in f_list:
                    fiber_bundles[-1].append(
                        h5f[group_name][str(fb)][str(f)][:].astype(float))
    else:
        raise TypeError(ext + ' is not implemented yet')

    return fiber_bundles


def save(file_name, fiber_bundles, group_name='fiber_bundles/', mode='w'):
    """
    Save fiberbundles configurations to a text file oder hdf5 file

    Parameters
    ----------
    file_name : string
        path to file
    fiber_bundles : list( list( (n,4)-array_like ) )
    group_name : string, optional
        path inside hdf5-file for fiber bundles
    mode : string, optional
        file mode of open() or h5py.File()
    """
    _, ext = os.path.splitext(file_name)

    if ext in ['.dat', '.txt']:
        with open(file_name, mode) as file:
            for fb, fiber_bundle in enumerate(fiber_bundles):
                for fiber in fiber_bundle:
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
        with h5py.File(file_name, mode) as h5f:
            if group_name[-1] is not '/':
                group_name = group_name + '/'

            for fb_i, fb in enumerate(fiber_bundles):
                grp_fb = h5f.create_group(group_name + str(fb_i) + '/')
                for i, f in enumerate(fb):
                    grp_fb[str(i)] = f[:, :]
    else:
        raise TypeError(ext + ' is not implemented yet')
