import numpy as np
import multiprocessing as mp
import ctypes

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


def parallel_rofl_exec(i):

    if not shr_mask[i]:
        directionmap[i] = 0
        inclmap[i] = 0
        trelmap[i] = 0
        dirdevmap[i] = 0
        incldevmap[i] = 0
        treldevmap[i] = 0
        funcmap[i] = 0
        itermap[i] = 0
    else:
        data = np.empty((dim[2] * dim[3]))
        for j in range(dim[2] * dim[3]):
            data[j] = shr_data[i * dim[2] * dim[3] + j]
        data = np.reshape(data, (dim[2], dim[3]))

        params, params_conf, func, n_iter = rofl_fit(direction_map[i], 6, 6,
                                                     dir_offset_.value,
                                                     tilt_angle_.value, data,
                                                     gain_.value,
                                                     grad_mode_.value, False)
        directionmap[i] = params[0]
        inclmap[i] = params[1]
        trelmap[i] = params[2]
        dirdevmap[i] = params_conf[0]
        incldevmap[i] = params_conf[1]
        treldevmap[i] = params_conf[2]
        funcmap[i] = func
        itermap[i] = n_iter


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

    global dim, shr_data, shr_mask
    global direction_map, directionmap, inclmap, trelmap
    global dirdevmap, incldevmap, treldevmap, funcmap, itermap
    global dir_offset_, tilt_angle_, gain_, grad_mode_

    gain_ = mp.Value('d', gain)
    dir_offset_ = mp.Value('d', dir_offset)
    tilt_angle_ = mp.Value('d', tilt_angle)
    grad_mode_ = mp.Value('i', grad_mode)

    data = np.array(data, dtype=float, copy=False)

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

    shr_data = mp.Array(ctypes.c_double, data.size, lock=False)
    loc_data = np.frombuffer(shr_data)
    loc_data[:] = np.moveaxis(data, [0, 1, 2, 3], [2, 0, 1, 3]).flatten()

    shr_mask = mp.Array(ctypes.c_int, mask.size, lock=False)
    loc_mask = np.frombuffer(shr_mask, dtype=np.int32)
    loc_mask[:] = mask.flatten()

    direction_map = mp.Array('d', loc_mask.size, lock=False)
    directionmap = mp.Array('d', loc_mask.size, lock=False)
    inclmap = mp.Array('d', loc_mask.size, lock=False)
    trelmap = mp.Array('d', loc_mask.size, lock=False)
    dirdevmap = mp.Array('d', loc_mask.size, lock=False)
    incldevmap = mp.Array('d', loc_mask.size, lock=False)
    treldevmap = mp.Array('d', loc_mask.size, lock=False)
    funcmap = mp.Array('d', loc_mask.size, lock=False)
    itermap = mp.Array('d', loc_mask.size, lock=False)

    loc_direction_map = np.frombuffer(direction_map)
    loc_direction_map[:] = epa.direction(data[0, :, :, :]).flatten()

    dim = [data.shape[1], data.shape[2], data.shape[0], data.shape[3]]

    with mp.Pool(processes=num_threads) as pool:
        pool.map(parallel_rofl_exec, range(loc_mask.size))

    return np.frombuffer(directionmap).reshape(
        mask.shape), np.frombuffer(inclmap).reshape(
            mask.shape), np.frombuffer(trelmap).reshape(
                mask.shape), np.frombuffer(dirdevmap).reshape(
                    mask.shape), np.frombuffer(incldevmap).reshape(
                        mask.shape), np.frombuffer(treldevmap).reshape(
                            mask.shape), np.frombuffer(funcmap).reshape(
                                mask.shape), np.frombuffer(itermap).reshape(
                                    mask.shape)
