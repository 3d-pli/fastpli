import numpy as np
import h5py
from tqdm import tqdm


def label_to_txt(data, file_name, gray_level, replace=[]):
    print("writing txt")
    with open(file_name, 'w') as file:
        file.write(
            str(data.shape[0]) + " " + str(data.shape[1]) + " " +
            str(data.shape[2]) + " " + str(gray_level) + "\n")

        # order has to be z,y,x in txt file
        data = np.transpose(data, (2, 1, 0))
        for mat in tqdm(data):
            if replace:
                mat_cp = mat.copy()
                for (key, val) in replace:
                    mat[mat_cp == key] = val

            np.savetxt(file, mat, delimiter='\t', fmt='%i')
            file.write('\n')


def label_to_binary(data, file_name, gray_level, replace=[]):
    data = np.array(data, dtype=np.int8)
    print("writing binary header")
    with open(file_name + '.header.txt', 'w') as file:
        file.write('{:12d}'.format(data.shape[0]) +
                   '{:12d}'.format(data.shape[1]) +
                   '{:12d}'.format(data.shape[2]) +
                   '{:12d}'.format(gray_level) +
                   ", ImageLx,ImageLy,ImageLz,ImageGray\n")
        file.write(file_name + ".binary , name of file in binary format\n")

    print("writing binary")
    # order has to be z,y,x in txt file
    data = np.transpose(data, (2, 1, 0))

    with open(file_name + '.binary', 'wb') as file:
        for mat in tqdm(data):
            if replace:
                mat_cp = mat.copy()
                for (key, val) in replace:
                    mat[mat_cp == key] = val

            mat.tofile(file)


def h5_to_txt(file_in, dset_name, file_out, gray_level, replace=[]):
    # TODO: chunking in case of mpi h5
    with h5py.File(file_in, 'r') as h5f:
        print("loading data")
        data = h5f[dset_name][:]
        label_to_txt(data, file_out, gray_level, replace)


def h5_to_binary(file_in, dset_name, file_out, gray_level, replace=[]):
    # TODO: chunking in case of mpi h5
    with h5py.File(file_in, 'r') as h5f:
        print("loading data")
        data = h5f[dset_name][:]
        label_to_binary(data, file_out, gray_level, replace)


if __name__ == "__main__":
    data = np.random.randint(5, size=(250, 250, 250))
    label_to_txt(data, "test.txt", 4)
    label_to_binary(data, "test", 4)

    # h5_to_binary(
    #     "fibers.gzip.h5",
    #     "tissue",
    #     "fibers.txt",
    #     2,
    #     replace=[(1, 2), (2, 1), (3, 0), (4, 1), (5, 2), (6, 2)])
