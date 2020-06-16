# -*- coding: utf-8 -*-
"""
Methods for manipulation of fiber_bundle objects
"""

import numpy as np
import copy

from . import fiber


def Rescale(fiber_bundle, scale, mod='all'):
    """
    Rescales fiber_bundle

    Parameters
    ----------
    fiber_bundle : [(,4)-array, ...]
        list of fibers
    scale : float
        scale factor
    mod : str, optional
        'all', 'points' or 'radii' will be scaled

    Returns
    -------
    res : [(,4)-array, ...]
        scaled fiber_bundle
    """

    fiber_bundle = copy.deepcopy(fiber_bundle)
    for i, _ in enumerate(fiber_bundle):
        fiber_bundle[i] = fiber.Rescale(fiber_bundle[i], scale, mod)
    return fiber_bundle


def Rotate(fiber_bundle, rot, offset=None):
    """
    Rotates fiber_bundle around offset

    Parameters
    ----------
    fiber_bundle : [(,4)-array, ...]
        list of fibers
    rot : (3,3)-array_like
        scale factor
    offset : 3d-array-array_like, optional
        offset for rotation center

    Returns
    -------
    res : [(,4)-array, ...]
        rotated fiber_bundle
    """

    fiber_bundle = copy.deepcopy(fiber_bundle)
    rot = np.array(rot, copy=False)
    if offset is not None:
        offset = np.array(offset, copy=False)
    for i, _ in enumerate(fiber_bundle):
        fiber_bundle[i] = fiber.Rotate(fiber_bundle[i], rot, offset)
    return fiber_bundle


def Translate(fiber_bundle, offset):
    """
    Translates fiber_bundle

    Parameters
    ----------
    fiber_bundle : [(,4)-array, ...]
        list of fibers
    offset : 3d-array-array_like
        offset to translate

    Returns
    -------
    res : [(,4)-array, ...]
        translated fiber_bundle
    """

    fiber_bundle = copy.deepcopy(fiber_bundle)
    offset = np.array(offset, copy=False)
    for i, _ in enumerate(fiber_bundle):
        fiber_bundle[i] = fiber.Translate(fiber_bundle[i], offset)
    return fiber_bundle


def ApplyFun(fiber_bundle, fun):
    """
    Applies function to fibers

    Parameters
    ----------
    fiber_bundle : [(,4)-array, ...]
        list of fibers
    fun : function

    Returns
    -------
    res : [(,4)-array, ...]
        translated fiber_bundle
    """

    fiber_bundle = copy.deepcopy(fiber_bundle)
    for i, _ in enumerate(fiber_bundle):
        fiber_bundle[i] = fun(fiber_bundle[i])
    return fiber_bundle


def ApplyFunToPosition(fiber_bundle, fun):
    """
    Applies function to fibers positions

    Parameters
    ----------
    fiber_bundle : [(,4)-array, ...]
        list of fibers
    fun : function

    Returns
    -------
    res : [(,4)-array, ...]
        translated fiber_bundle
    """

    fiber_bundle = copy.deepcopy(fiber_bundle)
    for i, _ in enumerate(fiber_bundle):
        fiber_bundle[i][:, :-1] = fun(fiber_bundle[i][:, :-1])
    return fiber_bundle


def ApplyFunToRadii(fiber_bundle, fun):
    """
    Applies function to fibers radii

    Parameters
    ----------
    fiber_bundle : [(,4)-array, ...]
        list of fibers
    fun : function

    Returns
    -------
    res : [(,4)-array, ...]
        translated fiber_bundle
    """

    fiber_bundle = copy.deepcopy(fiber_bundle)
    for i, _ in enumerate(fiber_bundle):
        fiber_bundle[i][:, -1] = fun(fiber_bundle[i][:, -1])
    return fiber_bundle


def Cut(fiber_bundle, voi):
    """
    Cut fiber into voi. The cutting process can create multiple fibers.
    It checks every cone_aabb if it overlapps with the voi.

    Parameters
    ----------
    fiber_bundle : [(,4)-array, ...]
        list of fibers
    voi : [xmin, ymin, zmin],[xmax,ymax,zmax]
        Volume of interest of which fibers to include. E.g. same as in
        Simulation

    Returns
    -------
    res : [(,4)-array, ...]
        cutted fiber_bundle
    """

    new_fiber_bundle = []
    for f in fiber_bundle:
        new_fiber_bundle.extend(fiber.Cut(f, voi))

    return new_fiber_bundle


def CutSphere(fiber_bundle, radius, center=[0, 0, 0]):
    """
    Cut fiber into sphere. The cutting process can create multiple fibers.
    It checks every cone_aabb if it overlapps with the sphere.

    Parameters
    ----------
    fiber_bundle : [(,4)-array, ...]
        list of fibers
    radius : float
        radius of cutting sphere
    center : 3d-array
        center of cutting sphere

    Returns
    -------
    res : [(,4)-array, ...]
        cutted fiber_bundle
    """

    new_fiber_bundle = []
    for f in fiber_bundle:
        new_fiber_bundle.extend(fiber.CutSphere(f, radius, center))

    return new_fiber_bundle
