import os
import time
import h5py
import numpy as np

from mpi4py import MPI
from fastpli.simulation import Simpli
import time
import h5py
import numpy as np
import sys
import os

from fastpli.simulation import Simpli

FILE_PATH = os.path.dirname(os.path.abspath(__file__))


# Read fiber data and prepair for PliGenerator
# TODO: write json -> parameter function

global_start = time.time()

simpli = Simpli()
# PliGeneration ###
simpli.pixel_size = 1
simpli.dim = [100, 100, 100]
simpli.ReadFiberFile('example/cube.h5')
simpli.SetFiberProperties([[(0.333, 0.004, 10, 'p'), (
    0.666, -0.004, 5, 'b'), (1.0, 0.004, 1, 'r')]])

# manipulation of fibers
# simpli.RotateVolumeAroundPoint(np.deg2rad(
#     20), np.deg2rad(-10), np.deg2rad(5), [10, -5, 7.5])
simpli.TranslateVolume([50, 50, 50])

with h5py.File(os.path.join(FILE_PATH, 'output_' + str(MPI.COMM_WORLD.Get_size()) + '.h5'), 'w', driver='mpio', comm=MPI.COMM_WORLD) as h5f:

    start = time.time()
    label_field, _, tissue_properties = simpli.GenerateTissue(only_label=True)
    label_field, vec_field, tissue_properties = simpli.GenerateTissue()
    dim_local, dim_offset = simpli.DimData()

    print(label_field.shape)

    dset = h5f.create_dataset(
        'tissue', (simpli.dim[0], simpli.dim[1], simpli.dim[2]), dtype=np.uint16)
    dset[dim_offset[0]:dim_offset[0] + dim_local[0], dim_offset[1]:dim_offset[1]
         + dim_local[1], dim_offset[2]:dim_offset[2] + dim_local[2]] = label_field

    dset = h5f.create_dataset(
        'vectorfield', (simpli.dim[0], simpli.dim[1], simpli.dim[2], 3), dtype=np.uint16)
    dset[dim_offset[0]:dim_offset[0] + dim_local[0], dim_offset[1]:dim_offset[1]
         + dim_local[1], dim_offset[2]:dim_offset[2] + dim_local[2], :] = vec_field

    end = time.time()

    print("TissueGeneration:", end - start)

'''
    # PliSimulation ###
    simpli.setup.filter_rotations = np.deg2rad([0, 30, 60, 90, 120, 150])
    simpli.setup.light_intensity = 26000
    simpli.setup.resolution = 1
    simpli.setup.untilt_sensor = True
    simpli.setup.wavelength = 525

    start = time.time()
    print("run_simulation: 0")
    image = simpli.run_simulation(
        label_field, vec_field, tissue_properties, 0, 0)
    h5f['data/0'] = image

    print("run_simulation: 1")
    image = simpli.run_simulation(
        label_field,
        vec_field,
        tissue_properties,
        np.deg2rad(5.5),
        np.deg2rad(0))
    h5f['data/1'] = image

    print("run_simulation: 2")
    image = simpli.run_simulation(
        label_field,
        vec_field,
        tissue_properties,
        np.deg2rad(5.5),
        np.deg2rad(90))
    h5f['data/2'] = image

    print("run_simulation: 3")
    image = simpli.run_simulation(
        label_field,
        vec_field,
        tissue_properties,
        np.deg2rad(5.5),
        np.deg2rad(180))
    h5f['data/3'] = image

    print("run_simulation: 4")
    image = simpli.run_simulation(
        label_field,
        vec_field,
        tissue_properties,
        np.deg2rad(5.5),
        np.deg2rad(270))
    h5f['data/4'] = image

    end = time.time()
    print("run_simulation:", end - start)
'''

global_end = time.time()
print("GlobalRuntime:", global_end - global_start)
