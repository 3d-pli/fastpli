# -*- coding: utf-8 -*-
"""
Solver Class
"""

from .__solver import _Solver

from ... import io
from ... import objects
from ... import tools
from ... import version

import numpy as np
import warnings
import os


class Solver(_Solver):
    """
    Solver Class for solving collisions between fibers
    """

    __is_frozen = False

    def __setattr__(self, key, value):
        if self.__is_frozen and not hasattr(self, key):
            raise TypeError("%r is a frozen class" % self)
        object.__setattr__(self, key, value)

    def __freeze(self):
        self.__is_frozen = True

    def __init__(self):
        super().__init__()

        self._drag = 0
        self._obj_min_radius = 0
        self._obj_mean_length = 0
        self._col_voi = None
        self._omp_num_threads = 1
        self._step_num = 0
        self.__display = None

        super()._set_omp_num_threads(self._omp_num_threads)

        self.__freeze()

    def get_dict(self):
        """ Get all member variables which are properties """
        members = dict()
        for key, value in self.__dict__.items():
            if key == '_cells_populations' or key == '_fiber_bundles':
                continue
            if key.startswith("_") and not key.startswith(
                    "__") and not key.startswith("_Solver"):
                if isinstance(value, np.ndarray):
                    members[key[1:]] = value.tolist()
                else:
                    members[key[1:]] = value
        return members

    def set_dict(self, input):
        """ Set dictionary of variables to class members """
        for key, value in input.items():
            if key.startswith("_"):
                raise ValueError("member variable cant be set directly")

            if value is not None:
                setattr(self, key, value)
            else:
                warnings.warn("None value in dict detected")

    @property
    def fiber_bundles(self):
        """ get/set fiber_bundles [[(,4)-array]] """
        return super()._get_fiber_bundles()

    @fiber_bundles.setter
    def fiber_bundles(self, fbs):
        fbs = objects.fiber_bundles.Cast(fbs)
        super()._set_fiber_bundles(fbs)

    @property
    def drag(self):
        """ drag value applied in each step, optional """
        return super()._get_parameters()[0]

    @drag.setter
    def drag(self, value):
        self._drag = value
        super()._set_parameters(self._drag, self._obj_min_radius,
                                self._obj_mean_length)

    @property
    def step_num(self):
        """ get/set number of applied steps """
        return self._step_num

    @drag.setter
    def step_num(self, value):
        self._step_num = int(value)

    @property
    def obj_min_radius(self):
        """ get/set minimal circular radius allowed for fiber """
        return super()._get_parameters()[1]

    @obj_min_radius.setter
    def obj_min_radius(self, value):
        self._obj_min_radius = value
        super()._set_parameters(self._drag, self._obj_min_radius,
                                self._obj_mean_length)

    @property
    def obj_mean_length(self):
        """ get/set mean value allowed for fiber segment """
        return super()._get_parameters()[2]

    @obj_mean_length.setter
    def obj_mean_length(self, value):
        self._obj_mean_length = value
        super()._set_parameters(self._drag, self._obj_min_radius,
                                self._obj_mean_length)

    @property
    def col_voi(self):
        """ get/set volume on witch the collision algorithm is applied """
        return self._col_voi

    @col_voi.setter
    def col_voi(self, voi):
        # TODO:
        if not isinstance(voi, tuple):
            raise TypeError("col_voi := (min, max)")
        self._col_voi = voi
        super()._set_col_voi(self._col_voi[0], self._col_voi[1])

    @property
    def omp_num_threads(self):
        """ get/set number of omp threads """
        return self._omp_num_threads

    @omp_num_threads.setter
    def omp_num_threads(self, num):
        self._omp_num_threads = int(num)
        self._omp_num_threads = super()._set_omp_num_threads(
            self._omp_num_threads)

    def step(self):
        """
        Applies collision solving algorithm for one step
        
        Returns True if solved
        """
        self._step_num += 1
        return super().step()

    def draw_scene(self, display=True):
        """ Draws model configuration in if OpenGl window can be created. """
        if self.__display is None:
            import platform
            if platform.system() == "Darwin":
                self.__display = display
            else:
                if "DISPLAY" in os.environ:
                    if os.environ['DISPLAY']:
                        self.__display = display
                    else:
                        warnings.warn("test_opengl: DISPLAY variable empty")
                        self.__display = False
                else:
                    warnings.warn("test_opengl: no DISPLAY variable detected")
                    self.__display = False

        if self.__display:
            super().draw_scene()

    def apply_boundary_conditions(self, n_max=10):
        """ Applies boundary conditions for n_max steps without collision solving. """
        if not isinstance(n_max, int) or n_max <= 0:
            raise TypeError("only integer > 0 allowed")

        super().apply_boundary_conditions(n_max)

        return self.fiber_bundles

    def save_parameter_h5(self, h5f, script=None):
        """ Saves class members without fiber_bundles in hdf5 file. """
        h5f.attrs['fastpli/version'] = version.__version__
        h5f.attrs['fastpli/compiler'] = version.__compiler__
        h5f.attrs['fastpli/libraries'] = version.__libraries__
        h5f.attrs['fastpli/pip_freeze'] = tools.helper.pip_freeze()
        h5f.attrs['fastpli/solver'] = str(self.get_dict())
        if script:
            h5f.attrs['script'] = script

    def save_h5(self, h5f, script=None):
        """ Saves class members in hdf5 file. """
        io.fiber_bundles.save_h5(h5f, self.fiber_bundles)
        self.save_parameter_h5(h5f, script)

    def load_h5(self, h5f):
        """ Loads class members from hdf5 file. """
        self.fiber_bundles = io.fiber_bundles.load_h5(h5f)
        self.set_dict(dict(eval(str(h5f.attrs['fastpli/solver']))))

        if h5f.attrs['fastpli/version'] != version.__version__:
            warnings.warn("__version__ changed")

        if h5f.attrs['fastpli/pip_freeze'] != tools.helper.pip_freeze():
            warnings.warn("pip_freeze changed")
