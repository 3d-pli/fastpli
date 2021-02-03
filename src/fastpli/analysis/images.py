# -*- coding: utf-8 -*-
"""
Methods for Fiber Orientation Maps and vectors
"""

# TODO: split into fom.py and vector.py

import numpy as np
import numba


@numba.njit()
def _vec_to_rgb(x, y, z):
    r = np.sqrt(x**2 + y**2 + z**2)
    x = int(x / r * 255)
    y = int(y / r * 255)
    z = int(z / r * 255)
    return np.array([np.abs(x), np.abs(y), np.abs(z)])


@numba.njit()
def _hsv_black_to_rgb_space(h, s, v):
    # images have to be saved in rgb space

    h = (h + 360) % 360

    hi = np.floor(h / 60)
    f = h / 60.0 - hi

    p = v * (1 - s)
    q = v * (1 - s * f)
    t = v * (1 - s * (1 - f))

    if hi == 1:
        r, g, b = q, v, p
    elif hi == 2:
        r, g, b = p, v, t
    elif hi == 3:
        r, g, b = p, q, v
    elif hi == 4:
        r, g, b = t, p, v
    elif hi == 5:
        r, g, b = v, p, q
    else:
        r, g, b = v, t, p

    return np.array((r * 255, g * 255, b * 255), np.int64)


@numba.njit()
def _orientation_to_hsv(directionValue, inclinationValue):
    h = 360.0 * np.abs(directionValue) / np.pi
    s = 1.0
    v = 1.0 - (2 * np.abs(inclinationValue) / np.pi)

    return _hsv_black_to_rgb_space(h, s, v)


@numba.njit()
def hsv_black_sphere(n=128):
    """
    Creates a hsv_black color sphere legend.

    Parameters
    ----------
    n : int
        length and height of resulting image

    Returns
    -------
    res : rgb-array
        resulting sphere
    """
    sphere = np.zeros((n, n, 3), dtype=np.uint8)

    for x in range(n):
        xx = (-n // 2 + x) / (n // 2)
        for y in range(n):
            yy = (-n // 2 + y) / (n // 2)
            if xx**2 + yy**2 > 1:
                continue
            zz = np.sqrt(1 - xx**2 - yy**2)

            direction = np.arctan2(yy, xx)
            inclination = np.arcsin(zz)

            if direction < 0:
                direction += np.pi

            sphere[x, y, :] = _orientation_to_hsv(direction, inclination)

    return sphere


@numba.njit()
def rgb_sphere(n=128):
    """
    Creates a rgb color sphere legend.

    Parameters
    ----------
    n : int
        length and height of resulting image

    Returns
    -------
    res : rgb-array
        resulting sphere
    """
    sphere = np.zeros((n, n, 3), dtype=np.uint8)

    for x in range(n):
        xx = (n // 2 - x) / (n // 2)
        for y in range(n):
            yy = (n // 2 - y) / (n // 2)
            if xx**2 + yy**2 > 1:
                continue
            zz = np.sqrt(1 - xx**2 - yy**2)

            sphere[x, y, :] = _vec_to_rgb(xx, yy, zz)

    return sphere


def unit_vectors(direction, inclination, mask=None):
    """
    Calculates the unit vector from direction and inclination

    Parameters
    ----------
    direction : 2d-array
        direction in radian
    inclination : 2d-array
        inclination in radian
    mask : 2d-array(bool), optional
        mask

    Returns
    -------
    res : 2d-array, 2d-array, 2d-array
        x-, y- and z-vector component in arrays
    """
    UnitX = np.sin(0.5 * np.pi - inclination) * np.cos(direction)
    UnitY = np.sin(0.5 * np.pi - inclination) * np.sin(direction)
    UnitZ = np.cos(0.5 * np.pi - inclination)

    if mask:
        _mask = np.logical_not(mask)
        UnitX[_mask] = 0
        UnitY[_mask] = 0
        UnitZ[_mask] = 0

    return UnitX, UnitY, UnitZ


def fom_hsv_black(direction, inclination, mask=None):
    """
    Calculates the fiber orientation map in hsv_black from direction and
    inclination

    Parameters
    ----------
    direction : 2d-array
        direction in radian
    inclination : 2d-array
        inclination in radian
    mask : 2d-array(bool), optional
        mask

    Returns
    -------
    res : rgb-array
        resulting hsv-black fom
    """
    if mask is None:
        mask = np.ones_like(direction, dtype=np.bool)

    hsv = np.zeros((mask.shape[0], mask.shape[1], 3), np.uint8)
    for x in range(mask.shape[0]):
        for y in range(mask.shape[1]):
            if not mask[x, y]:
                continue

            hsv[x, y, :] = _orientation_to_hsv(direction[x, y], inclination[x,
                                                                            y])
    return hsv


def fom_rgb(direction, inclination, mask=None):
    """
    Calculates the fiber orientation map in hsv_black from direction and
    inclination

    Parameters
    ----------
    direction : 2d-array
        direction in radian
    inclination : 2d-array
        inclination in radian
    mask : 2d-array(bool), optional
        mask

    Returns
    -------
    res : rgb-array
        resulting hsv-black fom
    """
    if mask is None:
        mask = np.ones_like(direction, dtype=np.bool)

    rgb = np.zeros((mask.shape[0], mask.shape[1], 3), np.uint8)
    for x in range(mask.shape[0]):
        for y in range(mask.shape[1]):
            if not mask[x, y]:
                continue

            rgb[x,
                y, :] = _vec_to_rgb(
                    np.sin(0.5 * np.pi - inclination[x, y]) *
                    np.cos(direction[x, y]),
                    np.sin(0.5 * np.pi - inclination[x, y]) *
                    np.sin(direction[x, y]),
                    np.cos(0.5 * np.pi - inclination[x, y]))
    return rgb
