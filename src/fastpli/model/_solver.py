from ._solver_cpp import _SolverCPP
from ..objects import Fiber

import numpy as np


class Solver(_SolverCPP):

    _drag = 0
    _obj_min_radius = 0
    _obj_mean_length = 0

    @property
    def fiber_bundles(self):
        fiber_bundles = super().get_fiber_bundles()
        for i in range(len(fiber_bundles)):
            for j in range(len(fiber_bundles[i])):
                fiber_bundles[i][j] = Fiber(
                    fiber_bundles[i][j].points,
                    fiber_bundles[i][j].radii)
        return fiber_bundles

    @fiber_bundles.setter
    def fiber_bundles(self, fbs):
        if not isinstance(fbs, list):
            raise TypeError("fbs is not a list")

        for fb_i, fb in enumerate(fbs):
            if not isinstance(fb, list):
                raise TypeError("fb is not a list")

            for f_i, f in enumerate(fb):
                if isinstance(f, Fiber):
                    continue
                elif isinstance(f, list):
                    f = np.array(f)

                if isinstance(f, np.ndarray):
                    if len(f.shape) is not 2 or f.shape[1] is not 4:
                        raise TypeError("fiber elements has to be of dim nx4")
                    fbs[fb_i][f_i] = Fiber(f[:, 0:-1], f[:, -1])
                else:
                    raise TypeError(
                        "fiber hast to be a objects.Fiber, 4d-list or 4d-array")

        # print(fbs)
        super().set_fiber_bundles(fbs)

    @property
    def drag(self):
        return super().get_parameters()[0]

    @drag.setter
    def drag(self, value):
        self._drag = value
        super(
        ).set_parameters(
            self._drag,
            self._obj_min_radius,
         self._obj_mean_length)

    @property
    def obj_min_radius(self):
        return super().get_parameters()[1]

    @obj_min_radius.setter
    def obj_min_radius(self, value):
        self._obj_min_radius = value
        super(
        ).set_parameters(
            self._drag,
            self._obj_min_radius,
         self._obj_mean_length)

    @property
    def obj_mean_length(self):
        return super().get_parameters()[2]

    @obj_mean_length.setter
    def obj_mean_length(self, value):
        self._obj_mean_length = value
        super(
        ).set_parameters(
            self._drag,
            self._obj_min_radius,
         self._obj_mean_length)

    @property
    def parameters(self):
        self._drag = super().get_parameters()[0]
        self._obj_min_radius = super().get_parameters()[1]
        self._obj_mean_length = super().get_parameters()[2]
        return super().get_parameters()

    @parameters.setter
    def parameters(self, tuple):
        super().set_parameters(tuple[0], tuple[1], tuple[2])
