import numpy as np


class _MPI:

    __is_frozen = False

    def __setattr__(self, key, value):
        if self.__is_frozen and not hasattr(self, key):
            raise TypeError(f'{self.__class__.__name__} is a frozen class')
        object.__setattr__(self, key, value)

    def __freeze(self):
        self.__is_frozen = True

    def __init__(self, mpi_comm):
        self._mpi_comm = mpi_comm
        self._rank = mpi_comm.Get_rank()
        self._num_p = mpi_comm.Get_size()

        self._gen_dim_local = None
        self._gen_dim_global = None
        self._gen_dim_offset = None
        self._optic_local_dim = None
        self._optic_global_dim = None

    def _set_gen_dim(self, loc, glob, offset):
        self._gen_dim_local = loc
        self._gen_dim_global = glob
        self._gen_dim_offset = offset

    def split_optic(self, data):
        self._optic_global_dim = data.shape
        self._optic_local_dim = np.array_split(np.arange(data.shape[0]),
                                               self._num_p)[self._rank]
        return np.array_split(data, self._num_p, axis=0)[self._rank]

    def save_split_h5(self, h5f, data_name, data):
        data = np.array(data, copy=False)
        global_dim = self._optic_global_dim[:data.ndim]
        dset = h5f.create_dataset(data_name, global_dim, dtype=data.dtype)
        if data.size:
            dset[self._optic_local_dim, :] = data

    def save_image_h5(self, h5f, data_name, data):
        data = data[self._gen_dim_offset[0]:self._gen_dim_offset[0] +
                    self._gen_dim_local[0],
                    self._gen_dim_offset[1]:self._gen_dim_offset[1] +
                    self._gen_dim_local[1], :]
        self.save_h5(h5f, data_name, data, 2)

    def save_h5(self, h5f, data_name, data, _lock_dim=None):
        """
        simpli can be seperated into different mpi processes.
        This function provides a parallel hdf5 io to save data
        inside the same h5-file.
        """

        data = np.array(data, copy=False)

        dset_dim = np.copy(self._gen_dim_global)
        if data.ndim < len(dset_dim):
            dset_dim = dset_dim[:len(data.shape)]
        if data.ndim > len(dset_dim):
            dset_dim = np.append(dset_dim, data.shape[3:])

        if _lock_dim:
            if isinstance(_lock_dim, int):
                _lock_dim = [_lock_dim]

            _lock_dim = list(_lock_dim)
            for i in _lock_dim:
                dset_dim[i] = data.shape[i]

        dset = h5f.create_dataset(data_name, dset_dim, dtype=data.dtype)

        if data.ndim == 2:
            if data.size * data.itemsize > 2 * (2**10)**3:  # 2 GB
                for i in range(data.shape[0]):
                    dset[i + self._gen_dim_offset[0],
                         self._gen_dim_offset[1]:self._gen_dim_offset[1] +
                         self._gen_dim_local[1]] = data[i, :]
            else:
                dset[self._gen_dim_offset[0]:self._gen_dim_offset[0] +
                     self._gen_dim_local[0],
                     self._gen_dim_offset[1]:self._gen_dim_offset[1] +
                     self._gen_dim_local[1]] = data

        elif data.ndim == 3:
            if data.size * data.itemsize > 2 * (2**10)**3:  # 2 GB
                for i in range(data.shape[0]):
                    dset[i + self._gen_dim_offset[0],
                         self._gen_dim_offset[1]:self._gen_dim_offset[1] +
                         self._gen_dim_local[1],
                         self._gen_dim_offset[2]:self._gen_dim_offset[2] +
                         self._gen_dim_local[2]] = data[i, :]
            else:
                dset[self._gen_dim_offset[0]:self._gen_dim_offset[0] +
                     self._gen_dim_local[0],
                     self._gen_dim_offset[1]:self._gen_dim_offset[1] +
                     self._gen_dim_local[1],
                     self._gen_dim_offset[2]:self._gen_dim_offset[2] +
                     self._gen_dim_local[2]] = data

        elif data.ndim > 3:
            if data.size * data.itemsize > 2 * (2**10)**3:  # 2 GB
                for i in range(data.shape[0]):
                    dset[i + self._gen_dim_offset[0],
                         self._gen_dim_offset[1]:self._gen_dim_offset[1] +
                         self._gen_dim_local[1],
                         self._gen_dim_offset[2]:self._gen_dim_offset[2] +
                         self._gen_dim_local[2], :] = data[i, :]
            else:
                dset[self._gen_dim_offset[0]:self._gen_dim_offset[0] +
                     self._gen_dim_local[0],
                     self._gen_dim_offset[1]:self._gen_dim_offset[1] +
                     self._gen_dim_local[1],
                     self._gen_dim_offset[2]:self._gen_dim_offset[2] +
                     self._gen_dim_local[2], :] = data

        else:
            raise TypeError('no compatible save_mpi_array_as_h5: ' + data_name)
