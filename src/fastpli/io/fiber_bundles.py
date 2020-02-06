import os
import numpy as np
import h5py

from .. import objects


def load(file_name, group_name='/'):
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

    if ext in ['.dat', '.txt']:
        with open(file_name, 'r') as f:
            fiber_bundles = load_dat(f)
    elif ext == '.h5':
        with h5py.File(file_name, 'r') as h5f:
            fiber_bundles = load_h5(h5f[group_name])
    else:
        raise TypeError(ext + ' is not implemented yet')

    return fiber_bundles


def load_dat(file):
    """
    Load fiberbundles configurations from a text file oder hdf5 file

    Parameters
    ----------
    file : file object
        opened file: with open(file_name, 'r') as file:

    Returns
    -------
    res : list(list(fiber)), fibers are (n,4)-arrays with (x,y,z,radii) for each fiber point
    """

    fiber = []
    fiber_bundles = [[]]
    for line in file:
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

    return fiber_bundles


def load_h5(h5f):
    """
    Load fiberbundles configurations from a hdf5 class

    Parameters
    ----------
    h5f : hdf5 class
        h5-file or group object

    Returns
    -------
    res : list(list(fiber)), fibers are (n,4)-arrays with (x,y,z,radii) for each fiber point
    """

    fiber_bundles = []
    fb_list = list(map(int, list(h5f.keys())))
    fb_list.sort()
    for fb in fb_list:
        fiber_bundles.append([])
        f_list = list(map(int, list(h5f[str(fb)].keys())))
        f_list.sort()
        for f in f_list:
            fiber_bundles[-1].append(h5f[str(fb)][str(f)][:].astype(float))

    return fiber_bundles


def save(file_name, fiber_bundles, group_name='/', mode='w-'):
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
        mode = 'a' if mode == 'w-' else mode
        with open(file_name, mode) as file:
            save_dat(file, fiber_bundles)
    elif ext == '.h5':
        mode = 'w-' if mode == 'a' else mode
        with h5py.File(file_name, mode) as h5f:
            save_h5(h5f.create_group(group_name), fiber_bundles)
    else:
        raise TypeError(ext + ' is not implemented yet')


def save_dat(file, fiber_bundles):
    """
    Save fiberbundles configurations to a text file oder hdf5 file

    Parameters
    ----------
    file_obj : file object
        opened file: with open(file_name, 'w-') as file:
    fiber_bundles : list( list( (n,4)-array_like ) )
    """

    fiber_bundles = objects.fiber_bundles.Cast(fiber_bundles)

    if not fiber_bundles:
        return

    for fb, fiber_bundle in enumerate(fiber_bundles):
        for fiber in fiber_bundle:
            if fiber.shape[1] != 4 or len(fiber.shape) != 2:
                raise TypeError('Wrong shape:', fiber.shape)
            for line in fiber:
                file.write(
                    str(line[0]) + " " + str(line[1]) + " " + str(line[2]) +
                    " " + str(line[3]) + "\n")
            file.write("\n")
        if fb != len(fiber_bundles) - 1:
            file.write("\n")


def save_h5(h5f, fiber_bundles):
    """
    Save fiberbundles configurations inside a hdf5 file

    Parameters
    ----------
    h5f : hdf5 class
        h5-file or group object
    fiber_bundles : list( list( (n,4)-array_like ) )
    """

    fiber_bundles = objects.fiber_bundles.Cast(fiber_bundles)
    if not fiber_bundles:
        return

    for fb_i, fb in enumerate(fiber_bundles):
        grp_fb = h5f.create_group(str(fb_i))
        for i, f in enumerate(fb):
            grp_fb[str(i)] = f[:, :]
