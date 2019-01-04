from ._helper.data_container import *

from numpy import prod as _prod

# TODO: add offset for mpi


def save_2_h5(dset, data, dim):

    if _prod(dim) != data.size():
        raise ValueError("_prod(dim) is not data.size")

    if data.size() > 2**32:
        n_chunk = _prod(dim[1:])
        data.set_data_chunk(n_chunk)
        for i in range(dim[0]):
            dset[i, :] = data.next_data_chunk(n_chunk).reshape(dim[1:])

    else:
        dset[:] = data.asarray().reshape(dim)
