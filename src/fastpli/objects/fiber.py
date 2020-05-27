# -*- coding: utf-8 -*-
"""
Methods for manipulation of fiber objects
"""

import numpy as np
import numba


def Rescale(fiber, scale, mod='all'):
    """
    Rescales fiber

    Parameters
    ----------
    fiber : (,4)-array
        fiber
    scale : float
        scale factor
    mod : str, optional
        'all', 'points' or 'radii' will be scaled

    Returns
    -------
    res : (,4)-array
        scaled fiber
    """

    fiber = np.array(fiber, copy=True)
    if mod == 'all':
        fiber *= scale
    elif mod == 'points':
        fiber[:, :3] *= scale
    elif mod == 'radii':
        fiber[:, -1] *= scale
    else:
        raise ValueError('mod = [all, points, radii]')
    return fiber


def Rotate(fiber, rot, offset=None):
    """
    Rotates fiber around offset

    Parameters
    ----------
    fiber : (,4)-array
        fiber
    rot : (3,3)-array_like
        scale factor
    offset : 3d-array-array_like, optional
        offset for rotation center

    Returns
    -------
    res : (,4)-array
        rotated fiber
    """

    rot = np.array(rot, copy=False)
    fiber = np.array(fiber, copy=True)
    if offset is None:
        fiber[:, :3] = np.dot(rot, fiber[:, :3].T).T
    else:
        offset = np.array(offset, copy=False)
        fiber[:, :3] = np.dot(rot, (fiber[:, :3] - offset).T).T + offset
    return fiber


def Translate(fiber, offset):
    """
    Translates fiber

    Parameters
    ----------
    fiber : (,4)-array
        fiber
    offset : 3d-array-array_like
        offset to translate

    Returns
    -------
    res : (,4)-array
        translated fiber
    """

    fiber = np.array(fiber, copy=True)
    offset = np.array(offset, copy=False)
    fiber[:, :3] += offset
    return fiber


@numba.njit(cache=True)
def _cone_aabb_in_aabb(c0, c1, vmin, vmax):
    c_min = np.array([
        min(c0[0] - c0[-1], c1[0] - c1[-1]),
        min(c0[1] - c0[-1], c1[1] - c1[-1]),
        min(c0[2] - c0[-1], c1[2] - c1[-1])
    ])

    c_max = np.array([
        min(c0[0] + c0[-1], c1[0] + c1[-1]),
        min(c0[1] + c0[-1], c1[1] + c1[-1]),
        min(c0[2] + c0[-1], c1[2] + c1[-1])
    ])

    for i in range(3):
        if c_min[i] > vmax[i] or c_max[i] < vmin[i]:
            return False
    return True


def Cut(fiber, voi):
    """
    Cut fiber into voi. The cutting process can create multiple fibers.
    It checks every cone_aabb if it overlapps with the voi.

    Parameters
    ----------
    fiber : (,4)-array
        fiber
    voi : [xmin, ymin, zmin],[xmax,ymax,zmax]
        Volume of interest of which fibers to include. E.g. same as in
        Simulation

    Returns
    -------
    res : [(,4)-array]
        cut fiber(s)
    """

    fibers = []
    fiber = np.array(fiber, copy=False)
    if fiber.ndim != 2:
        raise (TypeError, "False fiber shape")

    start = 0
    voi = np.array(voi)
    for i in range(fiber.shape[0] - 1):
        if not _cone_aabb_in_aabb(fiber[i, :], fiber[i + 1, :], voi[0], voi[1]):
            if start != i:
                fibers.append(fiber[start:i + 1])
            start = i + 1

    if start != i + 1:
        fibers.append(fiber[start:])

    return fibers
