# -*- coding: utf-8 -*-
"""
Simpli Class
"""

import warnings

from ._simpli import _Simpli


def Simpli(gpu=False, mpi_comm=None):
    if gpu:
        if mpi_comm:
            warnings.warn("no mpi communication implemented yet", UserWarning)
        pass
    else:
        return _Simpli(mpi_comm)
