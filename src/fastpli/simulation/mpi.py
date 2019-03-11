#!/usr/bin/python3

import numpy as np
from mpi4py import MPI

# TODO: integrade h5 get chunk for data > 2gb


def get_local(dim, fiber_bundles, max_theta):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    global_dim = dim.copy()

    # calculate dimension offset form rank
    local_dim_size, my_coords, _ = _local_dimensions(dim, rank, size)

    overlap = int(np.round(np.abs(dim[2] * np.tan(max_theta))))

    offset = local_dim_size * my_coords

    # add overlap offset
    offset -= overlap
    for fb in fiber_bundles:
        fb.translate_fiber_bundle([-offset[0], -offset[1], 0])

    dim[0] = local_dim_size[0]
    dim[1] = local_dim_size[1]

    if (offset[0] + overlap) + dim[0] > global_dim[0]:
        dim[0] = global_dim[0] - (offset[0] + overlap)
    if (offset[1] + overlap) + dim[1] > global_dim[1]:
        dim[1] = global_dim[1] - (offset[1] + overlap)

    dim[0] += 2 * overlap
    dim[1] += 2 * overlap

    return dim, fiber_bundles, offset, overlap


def crop_data(data, d):

    if data.ndim == 2:
        t_data = data[d:-d, d:-d]
    elif data.ndim == 3:
        t_data = data[d:-d, d:-d:]
    elif data.ndim == 4:
        t_data = data[d:-d, d:-d:, :]
    else:
        raise Exception("ERROR: not supported data type: " + str(data.ndim))

    return t_data


def _local_dimensions(dim, rank, size):

    g_dim = dim[0:2]
    min_area = float('Inf')

    for x in range(1, dim[0] + 1):
        for y in range(1, dim[1] + 1):
            if x * y != size:
                continue

            x = np.ceil(g_dim[0] / float(x))
            y = np.ceil(g_dim[1] / float(y))

            area = x + y

            if area < min_area:
                local_dim_size = np.array((x, y), np.int)
                min_area = area

    global_coords = np.ceil(g_dim / local_dim_size).astype(int)

    my_coords = np.array((rank % global_coords[0], rank // global_coords[0]),
                         np.int)

    return local_dim_size, my_coords, global_coords
