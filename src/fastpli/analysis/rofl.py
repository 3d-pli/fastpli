import numpy as np
import pymp

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

    direction = epa.direction(data)

    params, params_conf, func, n_iter = rofl_fit(direction, 6, 6, dir_offset,
                                                 tilt_angle, data, gain,
                                                 grad_mode, False)

    return params[0], params[1], params[2], params_conf[0], params_conf[
        1], params_conf[2], func, n_iter


def rofl_stack(data,
               tilt_angle=np.deg2rad(5.5),
               gain=3,
               dir_offset=0,
               mask=None,
               num_threads=2,
               grad_mode=False):
    '''
    data: np.array([tilts,x,y,stack])
    '''

    data = np.array(data, copy=False)

    if mask is None:
        mask = np.ones((data.shape[1], data.shape[2]), bool)

    if data.ndim != 4:
        raise TypeError("data: np.array([tilts,x,y,stack])")

    if data.shape[0] != 5:
        raise ValueError("data need 1 + 4 measurements")

    if data.shape[-1] <= 3:
        raise ValueError("data needs at least 3 equidistand rotations")

    if gain <= 0:
        raise ValueError("rofl gain <= 0")

    data_shr = pymp.shared.array(data.shape, data.dtype)
    data_shr = data

    direction_map = pymp.shared.array(data.shape[1:3], float)
    direction_map = pymp.shared.array(data.shape[1:3], float)
    direction_map[:, :] = epa.direction(data[0, :, :, :])

    directionmap = pymp.shared.array(data.shape[1:3], float)
    inclmap = pymp.shared.array(data.shape[1:3], float)
    trelmap = pymp.shared.array(data.shape[1:3], float)
    dirdevmap = pymp.shared.array(data.shape[1:3], float)
    incldevmap = pymp.shared.array(data.shape[1:3], float)
    treldevmap = pymp.shared.array(data.shape[1:3], float)
    funcmap = pymp.shared.array(data.shape[1:3], float)
    itermap = pymp.shared.array(data.shape[1:3], float)

    with pymp.Parallel(num_threads) as p:
        for x in p.range(data.shape[1]):
            for y in range(data.shape[2]):
                if not mask[x, y]:
                    directionmap[x, y] = 0
                    inclmap[x, y] = 0
                    trelmap[x, y] = 0
                    dirdevmap[x, y] = 0
                    incldevmap[x, y] = 0
                    treldevmap[x, y] = 0
                    funcmap[x, y] = 0
                    itermap[x, y] = 0
                else:
                    params, params_conf, func, n_iter = rofl_fit(
                        direction_map[x, y], 6, 6, dir_offset, tilt_angle,
                        data_shr[:, x, y, :], gain, grad_mode, False)
                    directionmap[x, y] = params[0]
                    inclmap[x, y] = params[1]
                    trelmap[x, y] = params[2]
                    dirdevmap[x, y] = params_conf[0]
                    incldevmap[x, y] = params_conf[1]
                    treldevmap[x, y] = params_conf[2]
                    funcmap[x, y] = func
                    itermap[x, y] = n_iter

    return np.array(directionmap), np.array(inclmap), np.array(
        trelmap), np.array(dirdevmap), np.array(incldevmap), np.array(
            treldevmap), np.array(funcmap), np.array(itermap)
