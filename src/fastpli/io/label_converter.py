import numpy as np
import h5py

try:
    from tqdm import tqdm
except ImportError as e:
    tqdm = lambda x: x


def label_to_txt(data, file_name, gray_level, replace=[]):
    print("writing txt")
    with open(file_name, 'w') as file:
        file.write(
            str(data.shape[0]) + " " + str(data.shape[1]) + " " +
            str(data.shape[2]) + " " + str(gray_level) + "\n")

        # order has to be z,y,x in txt file
        print("transposing data")
        data = np.ascontiguousarray(
            np.transpose(data.astype(np.int8), (2, 1, 0)))

        print("writing txt")
        for mat in tqdm(data):
            if replace:
                mat_cp = mat.copy()
                for (key, val) in replace:
                    mat[mat_cp == key] = val

            np.savetxt(file, mat, delimiter='\t', fmt='%i')
            file.write('\n')


def label_to_binary(data, file_name, gray_level, replace=[]):
    print("writing binary header")
    with open(file_name + '.header.txt', 'w') as file:
        file.write('{:12d}'.format(data.shape[0]) +
                   '{:12d}'.format(data.shape[1]) +
                   '{:12d}'.format(data.shape[2]) +
                   '{:12d}'.format(gray_level) +
                   ", ImageLx,ImageLy,ImageLz,ImageGray\n")
        file.write(file_name + ".binary , name of file in binary format\n")

    # order has to be z,y,x in txt file
    print("transposing data")
    data = np.ascontiguousarray(np.transpose(data.astype(np.int8), (2, 1, 0)))

    print("writing binary")
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
