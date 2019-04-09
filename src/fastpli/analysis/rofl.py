import numpy as np

from ._ROFL_with_jacobi import _execute_fit as rofl_fit
from . import epa


def rofl(data, tilt_angle=np.deg2rad(5.5), gain_value=3):
    '''
    data: np.array([tilts,stack])
    '''

    data = np.array(data)

    if len(data.shape) != 2:
        raise TypeError("data: np.array([tilts,stack])")

    if data.shape[0] != 5:
        raise ValueError("data need 1 + 4 measurements")

    if data.shape[1] <= 3:
        raise ValueError("data needs at least 3 equidistand rotations")

    direction = epa.direction(data)

    params, params_conf, func, n_iter = _execute_fit(
        direction, 6, 6, tilt_angle, data[:, x, y, :], gain_value)

    return params[0], params[1], params[2], params_conf[0], params_conf[
        1], params_conf[2], func, n_iter


def map(data, tilt_angle=np.deg2rad(5.5), gain_value=3, mask=None):
    '''
    data: np.array([tilts,x,y,stack])
    '''

    data = np.array(data)
    if mask is None:
        mask = np.ones((data.shape[1], data.shape[2]), bool)

    if len(data.shape) != 4:
        raise TypeError("data: np.array([tilts,x,y,stack])")

    if data.shape[0] != 5:
        raise ValueError("data need 1 + 4 measurements")

    if data.shape[-1] <= 3:
        raise ValueError("data needs at least 3 equidistand rotations")

    direction_map = epa.direction_map(data[0, :, :, :])

    directionmap = np.zeros_like(direction_map)
    inclmap = np.zeros_like(direction_map)
    trelmap = np.zeros_like(direction_map)
    dirdevmap = np.zeros_like(direction_map)
    incldevmap = np.zeros_like(direction_map)
    treldevmap = np.zeros_like(direction_map)
    funcmap = np.zeros_like(direction_map)
    itermap = np.zeros_like(direction_map)

    for x in range(data.shape[1]):
        for y in range(data.shape[2]):
            if not mask[x, y]:
                continue

            params, params_conf, func, n_iter = rofl_fit(
                direction_map[x, y], 6, 6, tilt_angle, data[:, x, y, :],
                gain_value)
            directionmap[x, y] = params[0]
            inclmap[x, y] = params[1]
            trelmap[x, y] = params[2]
            dirdevmap[x, y] = params_conf[0]
            incldevmap[x, y] = params_conf[1]
            treldevmap[x, y] = params_conf[2]
            funcmap[x, y] = func
            itermap[x, y] = n_iter

    return directionmap, inclmap, trelmap, dirdevmap, incldevmap, treldevmap, funcmap, itermap
