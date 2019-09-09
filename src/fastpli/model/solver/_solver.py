from .__solver import _Solver

from ...version import __version__

import numpy as np


class Solver(_Solver):

    __isfrozen = False

    def __setattr__(self, key, value):
        if self.__isfrozen and not hasattr(self, key):
            raise TypeError("%r is a frozen class" % self)
        object.__setattr__(self, key, value)

    def _freeze(self):
        self.__isfrozen = True

    def __init__(self):
        super().__init__()

        self._drag = 0
        self._obj_min_radius = 0
        self._obj_mean_length = 0
        self._col_voi = None
        self._omp_num_threads = 1

        super()._set_omp_num_threads(self._omp_num_threads)

        self._freeze()

    def as_dict(self):
        return {
            'version': __version__,
            'drag': self._drag,
            'obj_min_radius': self._obj_min_radius,
            'obj_mean_length': self._obj_mean_length,
            'col_voi': self._col_voi,
            'omp_num_threads': self._omp_num_threads
        }

    @property
    def fiber_bundles(self):
        return super()._get_fiber_bundles()

    @fiber_bundles.setter
    def fiber_bundles(self, fbs):
        if not isinstance(fbs, (list, tuple)):
            raise TypeError("fbs is not a list")

        for fb_i, fb in enumerate(fbs):
            if not isinstance(fb, (list, tuple)):
                raise TypeError("fb is not a list")

            for f_i, f in enumerate(fb):
                f = np.array(f, dtype=np.float64, copy=False)

                if len(f.shape) is not 2 or f.shape[1] is not 4:
                    raise TypeError("fiber elements has to be of dim nx4")

                fbs[fb_i][f_i] = f

        super()._set_fiber_bundles(fbs)

    @property
    def drag(self):
        return super()._get_parameters()[0]

    @drag.setter
    def drag(self, value):
        self._drag = value
        super()._set_parameters(self._drag, self._obj_min_radius,
                                self._obj_mean_length)

    @property
    def obj_min_radius(self):
        return super()._get_parameters()[1]

    @obj_min_radius.setter
    def obj_min_radius(self, value):
        self._obj_min_radius = value
        super()._set_parameters(self._drag, self._obj_min_radius,
                                self._obj_mean_length)

    @property
    def obj_mean_length(self):
        return super()._get_parameters()[2]

    @obj_mean_length.setter
    def obj_mean_length(self, value):
        self._obj_mean_length = value
        super()._set_parameters(self._drag, self._obj_min_radius,
                                self._obj_mean_length)

    @property
    def parameters(self):
        parameters = super()._get_parameters()
        self._drag = parameters[0]
        self._obj_min_radius = parameters[1]
        self._obj_mean_length = parameters[2]
        return parameters

    @property
    def col_voi(self):
        return self._col_voi

    @col_voi.setter
    def col_voi(self, voi):
        if not isinstance(voi, tuple):
            raise TypeError("col_voi := (min, max)")
        self._col_voi = voi
        super()._set_col_voi(self._col_voi[0], self._col_voi[1])

    @property
    def omp_num_threads(self):
        return self._omp_num_threads

    @omp_num_threads.setter
    def omp_num_threads(self, num):
        self._omp_num_threads = int(num)
        self._omp_num_threads = super()._set_omp_num_threads(
            self._omp_num_threads)
