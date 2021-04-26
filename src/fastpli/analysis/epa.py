# -*- coding: utf-8 -*-
"""
Analyse Methods for PLI signals
"""

import numpy as np


def epa(data):
    """
    Calculates modalities for a PLI image sequence

    Parameters
    ----------
    data : (x,y,rho)-array_like
        rho index must be equidistance between [0,180) degree

    Returns
    -------
    res : transmittance, direction, retardation
    """

    data = np.array(data, copy=False)
    n = data.shape[-1]

    dtype = np.float32 if data.itemsize <= 4 else np.float64
    rho_2 = 2 * np.linspace(0, np.pi, n, False, dtype=dtype)

    a0 = np.sum(data, -1) / n
    a1 = 2 * np.sum(data * np.sin(rho_2), -1) / n
    b1 = 2 * np.sum(data * np.cos(rho_2), -1) / n

    t = 2 * a0
    d = 0.5 * np.arctan2(-b1, a1) + np.pi
    r = np.sqrt(a1 * a1 + b1 * b1) / (a0 + 1e-16)

    d = d % np.pi

    # TODO: d = 0.5 * np.arctan2(a1, -b1) + np.pi without d = d % np.pi

    return t, d, r


def direction(data):
    """
    Calculates direction map for a PLI image sequence

    Parameters
    ----------
    data : array_like
        (x,y,rho)-array, rho index must be equidistance between [0,180) degree

    Returns
    -------
    res : direction
    """

    data = np.array(data, copy=False)
    n = data.shape[-1]

    dtype = np.float32 if data.itemsize <= 4 else np.float64
    rho_2 = 2 * np.linspace(0, np.pi, n, False, dtype=dtype)

    a1 = 2 * np.sum(data * np.sin(rho_2), -1) / n
    b1 = 2 * np.sum(data * np.cos(rho_2), -1) / n

    d = 0.5 * np.arctan2(-b1, a1) + np.pi

    d = d % np.pi

    return d
