from .__solver import __Solver
from ..objects import Fiber

import numpy as np


class Solver(__Solver):

    _drag = 0
    _obj_min_radius = 0
    _obj_mean_length = 0
    _col_voi = None

    @property
    def fiber_bundles(self):
        fiber_bundles = super()._get_fiber_bundles()
        for i in range(len(fiber_bundles)):
            for j in range(len(fiber_bundles[i])):
                fiber_bundles[i][j] = Fiber(
                    fiber_bundles[i][j].points,
                    fiber_bundles[i][j].radii)
        return fiber_bundles

    @fiber_bundles.setter
    def fiber_bundles(self, fbs):
        if not isinstance(fbs, (list, tuple)):
            raise TypeError("fbs is not a list")

        for fb_i, fb in enumerate(fbs):
            if not isinstance(fb, (list, tuple)):
                raise TypeError("fb is not a list")

            for f_i, f in enumerate(fb):
                if isinstance(f, Fiber):
                    continue
                elif isinstance(f, (list, tuple)):
                    f = np.array(f)

                if isinstance(f, np.ndarray):
                    if len(f.shape) is not 2 or f.shape[1] is not 4:
                        raise TypeError("fiber elements has to be of dim nx4")
                    fbs[fb_i][f_i] = Fiber(f[:, 0:-1], f[:, -1])
                else:
                    raise TypeError(
                        "fiber hast to be a objects.Fiber, 4d-list or 4d-array")

            super()._set_fiber_bundles(fbs)

    @property
    def drag(self):
        return super()._get_parameters()[0]

    @drag.setter
    def drag(self, value):
        self._drag = value
        super()._set_parameters(
            self._drag,
            self._obj_min_radius,
            self._obj_mean_length)

    @property
    def obj_min_radius(self):
        return super()._get_parameters()[1]

    @obj_min_radius.setter
    def obj_min_radius(self, value):
        self._obj_min_radius = value
        super()._set_parameters(
            self._drag,
            self._obj_min_radius,
            self._obj_mean_length)

    @property
    def obj_mean_length(self):
        return super()._get_parameters()[2]

    @obj_mean_length.setter
    def obj_mean_length(self, value):
        self._obj_mean_length = value
        super()._set_parameters(
            self._drag,
            self._obj_min_radius,
            self._obj_mean_length)

    @property
    def parameters(self):
        parameters = super()._get_parameters()
        self._drag = parameters[0]
        self._obj_min_radius = parameters[1]
        self._obj_mean_length = parameters[2]
        return parameters

    @parameters.setter
    def parameters(self, tuple):
        super()._set_parameters(tuple[0], tuple[1], tuple[2])

    def set_parameters(self, drag=0, obj_min_radius=0, obj_mean_length=0):
        super()._set_parameters(drag, obj_min_radius, obj_mean_length)

    @property
    def col_voi(self):
        return self._col_voi

    @col_voi.setter
    def col_voi(self, voi):
        if not isinstance(voi, tuple):
            raise TypeError("col_voi := (min, max)")
        self._col_voi = voi
        super()._set_col_voi(self._col_voi[0], self._col_voi[1])
