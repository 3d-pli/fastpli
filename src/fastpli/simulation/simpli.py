# -*- coding: utf-8 -*-
"""
Simpli Class
"""

import warnings

from . import simpli_cpu


def Simpli(gpu=False, mpi_comm=None):
    if gpu:
        if mpi_comm:
            warnings.warn('no mpi communication implemented yet', UserWarning)
        raise ValueError('GPU Version will be added shortly')
    else:
        return simpli_cpu.Simpli(mpi_comm)
