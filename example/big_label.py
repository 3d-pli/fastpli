import fastpli.io
import fastpli.model
import fastpli.simulation

import h5py
import numpy as np
import os
from tqdm import tqdm

print(
    "WARNING: chunking can change the label_field because of numerical precision"
)

INPUT_FILE_FIBER = '/tmp/test.dat'
OUTPUT_FILE_TISSUE = '/tmp/test_chunked.h5'

# tissue generation
simpli = fastpli.simulation.Simpli()
simpli.debug = False

simpli.pixel_size = 0.025  # in mu-meter
simpli.dim = [200, 200, 200]  # num of voxels
simpli.dim_origin = [25, 15, 15]  # in mu-meter
# or: simpli.voi = [25, 30, 15, 20, 15, 20] in mu-meter

simpli.fiber_bundles = fastpli.io.fiber.read(INPUT_FILE_FIBER)
simpli.fiber_bundles_properties = [[(0.6, 0, 0, 'p'), (0.8, 0, 0, 'p'),
                                    (1, 0, 0, 'p')]]

with h5py.File(OUTPUT_FILE_TISSUE, 'w') as h5f:
    dset = h5f.create_dataset('tissue',
                              simpli.dim,
                              dtype=np.uint16,
                              compression="gzip")

    # chunking voi in x-direction
    X_STEP = 1
    x_voi = (simpli.get_voi()[0], simpli.get_voi()[1])
    x_list = [x for x in np.arange(x_voi[0], x_voi[1], X_STEP)]
    simpli.dim = [simpli.dim[0] / len(x_list), simpli.dim[1], simpli.dim[2]]

    for x_id, x_pos in enumerate(tqdm(x_list)):
        simpli._dim_origin[0] = x_pos
        label_field, _, _ = simpli.GenerateTissue(only_label=True)
        x_offset = x_id * simpli.dim[0]
        dset[x_offset:x_offset + label_field.shape[0], :, :] = label_field
