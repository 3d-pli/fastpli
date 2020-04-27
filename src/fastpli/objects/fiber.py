# -*- coding: utf-8 -*-
"""
Methods for manipulation of fiber objects
"""

import numpy as np


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
