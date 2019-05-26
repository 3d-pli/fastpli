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
FILE_NAME = '/tmp/fastpli.example.simpli_mpi.' + str(
    MPI.COMM_WORLD.Get_size()) + '.h5'
with h5py.File(FILE_NAME, 'w', driver='mpio', comm=MPI.COMM_WORLD) as h5f:

    # Read fiber data and prepair for PliGenerator
    # TODO: write json -> parameter function

    global_start = time.time()

    simpli = Simpli(MPI.COMM_WORLD)

    # PliGeneration ###
    simpli.pixel_size = 1
    simpli.set_voi([-50, 50, -50, 50, -50, 50])
    simpli.ReadFiberFile(os.path.join(FILE_PATH, 'cube.h5'))
    simpli.fiber_bundles_properties = ([[(0.333, 0.004, 10, 'p'),
                                         (0.666, -0.004, 5, 'b'),
                                         (1.0, 0.004, 1, 'r')]])

    # manipulation of fibers
    # simpli.RotateVolumeAroundPoint(
    #     np.deg2rad(20), np.deg2rad(-10), np.deg2rad(5), [10, -5, 7.5])
    # simpli.TranslateVolume([50, 50, 50])

    start = time.time()
    label_field, vec_field, tissue_properties = simpli.GenerateTissue(
        only_label=False)
    dim_local, dim_offset = simpli.DimData()

    simpli.SaveAsH5(h5f, label_field, "tissue")
    simpli.SaveAsH5(h5f, vec_field, "vectorfield")

    end = time.time()
    print("TissueGeneration:", end - start)

    # PliSimulation ###
    simpli.filter_rotations = np.deg2rad([0, 30, 60, 90, 120, 150])
    simpli.light_intensity = 26000
    simpli.untilt_sensor = True
    simpli.wavelength = 525

    start = time.time()
    print("RunSimulation: 0")
    image = simpli.RunSimulation(label_field, vec_field, tissue_properties, 0,
                                 0)
    simpli.SaveAsH5(h5f, image, 'data/0')

    print("RunSimulation: 1")
    image = simpli.RunSimulation(label_field, vec_field, tissue_properties,
                                 np.deg2rad(5.5), np.deg2rad(0))
    simpli.SaveAsH5(h5f, image, 'data/1')

    print("RunSimulation: 2")
    image = simpli.RunSimulation(label_field, vec_field, tissue_properties,
                                 np.deg2rad(5.5), np.deg2rad(90))
    simpli.SaveAsH5(h5f, image, 'data/2')

    print("RunSimulation: 3")
    image = simpli.RunSimulation(label_field, vec_field, tissue_properties,
                                 np.deg2rad(5.5), np.deg2rad(180))
    simpli.SaveAsH5(h5f, image, 'data/3')

    print("RunSimulation: 4")
    image = simpli.RunSimulation(label_field, vec_field, tissue_properties,
                                 np.deg2rad(5.5), np.deg2rad(270))
    simpli.SaveAsH5(h5f, image, 'data/4')

    print("RunSimulation: 5")
    image = simpli.RunSimulation(label_field, vec_field, tissue_properties,
                                 np.deg2rad(5.5), np.deg2rad(42))
    simpli.SaveAsH5(h5f, image, 'data/5')

    end = time.time()
    print("RunSimulation:", end - start)

    global_end = time.time()
    print("GlobalRuntime:", global_end - global_start)
