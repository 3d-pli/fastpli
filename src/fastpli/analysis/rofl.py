import numpy as np

from ._ROFL_with_jacobi import _execute_fit as rofl_fit
from . import epa


def rofl(data,
         tilt_angle=np.deg2rad(5.5),
         gain=3,
         dir_offset=0,
         grad_mode=False):
    '''
    data: np.array([tilts,stack])
    '''

    data = np.array(data, copy=False)

    if data.ndim != 2:
        raise TypeError("data: np.array([tilts,stack])")

    if data.shape[0] != 5:
        raise ValueError("data need 1 + 4 measurements")

    if data.shape[1] <= 3:
        raise ValueError("data needs at least 3 equidistand rotations")

    if gain <= 0:
        raise ValueError("rofl gain <= 0")

    direction = epa.direction(data[0, :])

    params, params_conf, func, n_iter = rofl_fit(direction, 6, 6, dir_offset,
                                                 tilt_angle, data, gain,
                                                 grad_mode, False)

    return params[0], params[1], params[2], params_conf[0], params_conf[
        1], params_conf[2], func, n_iter
