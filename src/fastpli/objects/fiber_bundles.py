# -*- coding: utf-8 -*-
"""
Methods for manipulation of fiber_bundles objects
"""

import numpy as np
import copy

from . import fiber


def Cast(fiber_bundles):
    """
    Cast objects into fiber_bundle object

    Parameters
    ----------
    fiber_bundles : [[(,4)-array_like, ...]-like]-like
        list of fiber_bundle

    Returns
    -------
    res : [[(,4)-array, ...]]
        fiber_bundles
    """

    if not fiber_bundles:
        return fiber_bundles

    if not isinstance(fiber_bundles, (list, tuple)):
        raise TypeError("fiber_bundles is not a list")

    for fb_i, fb in enumerate(fiber_bundles):
        if not isinstance(fb, (list, tuple)):
            raise TypeError("fiber_bundle is not a list")

        for f_i, f in enumerate(fb):
            fiber_bundles[fb_i][f_i] = np.array(f, dtype=float)

            if fiber_bundles[fb_i][f_i].ndim != 2 or fiber_bundles[fb_i][
                    f_i].shape[1] != 4:
                raise TypeError("fiber elements has to be of shape nx4")

    return fiber_bundles


def Rescale(fiber_bundles, scale, mod='all'):
    """
    Rescales fiber_bundles

    Parameters
    ----------
    fiber_bundles : [[(,4)-array, ...]]
        list of fiber_bundle
    scale : float
        scale factor
    mod : str, optional
        'all', 'points' or 'radii' will be scaled

    Returns
    -------
    res : [[(,4)-array, ...]]
        scaled fiber_bundles
    """

    fiber_bundles = copy.deepcopy(fiber_bundles)
    for j, fb in enumerate(fiber_bundles):
        for i, _ in enumerate(fb):
            fiber_bundles[j][i] = fiber.Rescale(fiber_bundles[j][i], scale, mod)
    return fiber_bundles


def Rotate(fiber_bundles, rot, offset=None):
    """
    Rotates fiber_bundles around offset

    Parameters
    ----------
    fiber_bundles : [[(,4)-array, ...]]
        list of fibers
    rot : (3,3)-array_like
        scale factor
    offset : 3d-array-array_like, optional
        offset for rotation center

    Returns
    -------
    res : [[(,4)-array, ...]]
        rotated fiber_bundles
    """

    fiber_bundles = copy.deepcopy(fiber_bundles)
    rot = np.array(rot, copy=False)
    if offset is not None:
        offset = np.array(offset, copy=False)
    for j, fb in enumerate(fiber_bundles):
        for i, _ in enumerate(fb):
            fiber_bundles[j][i] = fiber.Rotate(fiber_bundles[j][i], rot, offset)
    return fiber_bundles


def Translate(fiber_bundles, offset):
    """
    Translates fiber_bundles

    Parameters
    ----------
    fiber_bundles : [[(,4)-array, ...]]
        list of fibers
    offset : 3d-array-array_like
        offset to translate

    Returns
    -------
    res : [[(,4)-array, ...]]
        translated fiber_bundles
    """

    fiber_bundles = copy.deepcopy(fiber_bundles)
    offset = np.array(offset, copy=False)
    for j, fb in enumerate(fiber_bundles):
        for i, _ in enumerate(fb):
            fiber_bundles[j][i] = fiber.Translate(fiber_bundles[j][i], offset)
    return fiber_bundles


def ApplyFun(fiber_bundles, fun):
    """
    Applies function to fibers

    Parameters
    ----------
    fiber_bundles : [[(,4)-array, ...]]
        list of fibers
    fun : function

    Returns
    -------
    res : [[(,4)-array, ...]]
        fiber_bundles
    """

    fiber_bundles = copy.deepcopy(fiber_bundles)
    for j, fb in enumerate(fiber_bundles):
        for i, _ in enumerate(fb):
            fiber_bundles[j][i] = fun(fiber_bundles[j][i])
    return fiber_bundles


def ApplyFunToPosition(fiber_bundles, fun):
    """
    Applies function to fibers positions

    Parameters
    ----------
    fiber_bundles : [[(,4)-array, ...]]
        list of fibers
    fun : function

    Returns
    -------
    res : [[(,4)-array, ...]]
        fiber_bundles
    """

    fiber_bundles = copy.deepcopy(fiber_bundles)
    for j, fb in enumerate(fiber_bundles):
        for i, _ in enumerate(fb):
            fiber_bundles[j][i][:, :-1] = fun(fiber_bundles[j][i][:, :-1])
    return fiber_bundles


def ApplyFunToRadii(fiber_bundles, fun):
    """
    Applies function to fibers radii

    Parameters
    ----------
    fiber_bundles : [[(,4)-array, ...]]
        list of fibers
    fun : function

    Returns
    -------
    res : [[(,4)-array, ...]]
        fiber_bundles
    """

    fiber_bundles = copy.deepcopy(fiber_bundles)
    for j, fb in enumerate(fiber_bundles):
        for i, _ in enumerate(fb):
            fiber_bundles[j][i][:, -1] = fun(fiber_bundles[j][i][:, -1])
    return fiber_bundles
