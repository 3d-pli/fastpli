import numpy as np
import h5py
from tqdm import tqdm

def label_to_txt(data, file_name, replace=[]):
    print(data.shape)

    print("writing txt")
    with open(file_name, 'w') as file:
        file.write( str(data.shape[0]) + " " + str(data.shape[1]) + " " + str(data.shape[2]) + " 2\n")

        # order has to be z,y,x in txt file
        data = np.transpose(data, (2,1,0))
        for mat in tqdm(data):
            if replace:
                mat_cp = mat.copy()
                for (key, val) in replace:
                    mat[mat_cp == key] = val

            np.savetxt(file, mat, delimiter='\t', fmt='%i')
            file.write('\n')

def h5_to_txt(file_in, dset_name, file_out, replace=[]):
    # TODO: chunking in case of mpi h5
    with h5py.File(file_in, 'r') as h5f:
        print("loading data")
        data = h5f[dset_name][:]
        label_to_txt(data, file_out, replace)
        # label_to_txt_old(data, file_out, replace)
