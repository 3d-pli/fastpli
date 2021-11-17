# -*- coding: utf-8 -*-
"""
Analyse Methods for tilted PLI signals with ROFL algorithm

Notes
-----
Algorithm public available at https://doi.org/10.3389/fnana.2018.00075

"""

import numpy as np

from ._ROFL_with_jacobi import _execute_fit as rofl_fit
from . import epa


def rofl(data,
         tilt_angle=np.deg2rad(5.5),
         gain=3,
         dir_offset=0,
         grad_mode=False):
    """
    Calculates modalities for a PLI image sequence

    Parameters
    ----------
    data : (tilts, rho)-array_like
        tilts must be [(0,0), (tilt_angle,0), (tilt_angle,pi/2),
        (tilt_angle,pi), (tilt_angle,3/2*pi)]
        rho must be equidistance between [0,pi)

    Returns
    -------
    res : direction, inclination, trel, ...
    """

    data = np.array(data, copy=False)
    dtype = np.float32 if data.itemsize <= 4 else np.float64
    data = data.astype(dtype)

    if data.ndim != 2:
        raise TypeError('data: np.array([tilts,stack])')

    if data.shape[0] != 5:
        raise ValueError('data need 1 + 4 measurements')

    if data.shape[1] <= 3:
        raise ValueError('data needs at least 3 equidistant rotations')

    if gain <= 1:
        raise ValueError('rofl gain <= 1')

    direction = epa.direction(data[0, :]) + dir_offset

    params, params_conf, func, n_iter = rofl_fit(direction, 6, 6, dir_offset,
                                                 tilt_angle, data, gain,
                                                 grad_mode, False)

    return params[0], params[1], params[2], params_conf[0], params_conf[
        1], params_conf[2], func, n_iter
