# -*- coding: utf-8 -*-
"""
fiber objects classes
"""

import copy

import numpy as np
import numba

from .layers import Layers
"""
--------------------------------------------------------------------------------
--------------------------------------FIBER-------------------------------------
--------------------------------------------------------------------------------
"""


def _convert_to_fiber(data, dtype):
    """ converts data into (n,4)-np.array"""

    if data is None:
        return np.empty((), dtype)
    if isinstance(data, Fiber):
        return data

    data = np.atleast_2d(np.array(data, dtype=dtype, copy=True))
    if data.ndim != 2:
        raise TypeError('Wrong shape: expected (n,4)')
    if data.shape[1] != 4:
        raise TypeError('Wrong shape: expected (n,4)')
    if not np.issubdtype(data.dtype, np.floating):
        raise TypeError('Wrong type: has to be floating')

    return data


class Fiber:
    __is_frozen = False

    def __setattr__(self, key, value):
        if self.__is_frozen and not hasattr(self, key):
            raise TypeError('%r is a frozen class' % self)
        object.__setattr__(self, key, value)

    def __freeze(self):
        self.__is_frozen = True

    def __init__(self, data, dtype=float):
        self._data = _convert_to_fiber(data, dtype)
        self.__freeze()

    def __getitem__(self, item):
        return self._data[item]

    def __setitem__(self, item, value):
        self._data[item] = value

    def __str__(self):
        return self._data.__str__()

    def __repr__(self):
        return self._data.__repr__()

    def copy(self):
        """ deep copy of class """
        return copy.deepcopy(self)

    def __iter__(self):
        return self._data.__iter__()

    def __next__(self):
        return self._data.__next__()

    def __len__(self):
        return self._data.shape[-1]

    @property
    def shape(self):
        """ returns shape of data """
        return self._data.shape

    @property
    def dtype(self):
        """ returns np.dtype of data """
        return self._data.dtype

    def as_array(self):
        """
        Returns copy data as np.array

        Returns
        -------
        res : (n,4)-array
            fiber as np.array
        """
        return self._data.copy()

    def cast(self, dtype):
        """
        Changes datatype to new type

        Parameters
        ----------
        dtype : type
            numpy types are fully supported

        Returns
        -------
        res : Fiber
            casted fiber
        """

        return Fiber(self, dtype)

    @property
    def points(self):
        """
        Returns xyz data as np.array

        Returns
        -------
        res : (n,3)-array
            fiber points as numpy array
        """
        return self._data[:, :-1]

    @points.setter
    def points(self, value):
        self._data[:, :-1] = value

    @property
    def radii(self):
        """
        Returns radii data as np.array

        Returns
        -------
        res : (n,1)-array
            fiber radii as numpy array
        """
        return self._data[:, -1]

    @radii.setter
    def radii(self, value):
        self._data[:, -1] = value

    def scale(self, scale, mode='all'):
        """
        Scales fiber

        Parameters
        ----------
        scale : float
            scale factor
        mode : str, optional
            'all', 'points' or 'radii' will be scaled

        Returns
        -------
        res : Fiber
            scaled Fiber
        """

        data = self._data.copy()

        if mode == 'all':
            data[:] *= scale
        elif mode == 'points':
            data[:, :3] *= scale
        elif mode == 'radii':
            data[:, -1] *= scale
        else:
            raise ValueError('mode = [all, points, radii]')

        return Fiber(data)

    def rotate(self, rot, offset=None):
        """
        Rotates fiber around offset

        Parameters
        ----------
        rot : (3,3)-array_like
            rotation matrix
        offset : 3d-array-array_like, optional
            offset for rotation center

        Returns
        -------
        res : Fiber
            rotated fiber
        """

        data = self._data.copy()

        rot = np.array(rot, copy=False)
        if offset is None:
            data[:, :3] = np.dot(rot, data[:, :3].T).T
        else:
            offset = np.array(offset, copy=False)
            data[:, :3] = np.dot(rot, (data[:, :3] - offset).T).T + offset

        return Fiber(data)

    def translate(self, offset):
        """
        Translates fiber

        Parameters
        ----------
        offset : 3d-array-array_like
            offset to translate

        Returns
        -------
        res : Fiber
            translated fiber
        """

        data = self._data.copy()

        offset = np.array(offset, copy=False)
        data[:, :3] += offset

        return Fiber(data)

    def apply(self, fun):
        """
        Applies function to fiber

        Parameters
        ----------
        fun : function

        Returns
        -------
        res : Fiber
            fun(fiber)
        """

        fiber = self.copy()
        fiber[:] = fun(fiber[:])
        return fiber

    def apply_to_points(self, fun):
        """
        Applies function to fiber positions

        Parameters
        ----------
        fun : function

        Returns
        -------
        res : Fiber
            fun(fiber)
        """

        fiber = self.copy()
        fiber[:, :-1] = fun(fiber[:, :-1])
        return fiber

    def apply_to_radii(self, fun):
        """
        Applies function to fiber radii

        Parameters
        ----------
        fun : function

        Returns
        -------
        res : Fiber
            fun(fiber)
        """

        fiber = self.copy()
        fiber[:, -1] = fun(fiber[:, -1])
        return fiber

    def cut(self, voi):
        """
        Cut fiber into voi. The cutting process can create multiple fibers.
        It checks every fiber_segment_aabb if it overlapps with the voi.

        Parameters
        ----------
        voi : [xmin, ymin, zmin],[xmax,ymax,zmax]
            Volume of interest of which fibers to include. E.g. same as in
            Simulation

        Returns
        -------
        res : FiberBundle
            cutted fiber into fiber_bundle
        """

        fibers = []

        start = 0
        voi = np.array(voi)
        for i in range(self._data.shape[0] - 1):
            if not _fiber_segment_aabb_in_aabb(
                    self._data[i, :], self._data[i + 1, :], voi[0], voi[1]):
                if start != i:
                    fibers.append(self._data[start:i + 1])
                start = i + 1

        if start != i + 1:
            fibers.append(self._data[start:])

        return FiberBundle(fibers)

    def cut_sphere(self, radius, center=(0, 0, 0)):
        """
        Cut fiber into sphere. The cutting process can create multiple fibers.
        It checks every fiber_segment_aabb if it overlapps with the sphere.

        Parameters
        ----------
        radius : float
            radius of cutting sphere
        center : 3d-array
            center of cutting sphere

        Returns
        -------
        res : FiberBundle
            cutted fiber_bundle
        """

        center = np.array(center, copy=False)

        fibers = []

        start = 0
        for i in range(self._data.shape[0] - 1):
            if not _fiber_segment_aabb_in_sphere(
                    self._data[i, :], self._data[i + 1, :], radius, center):
                if start != i:
                    fibers.append(self._data[start:i + 1])
                start = i + 1

        if start != i + 1:
            fibers.append(self._data[start:])

        return FiberBundle(fibers)


@numba.njit(cache=True)
def _fiber_segment_aabb_in_aabb(c0, c1, vmin, vmax):
    c_min = np.array([
        min(c0[0] - c0[-1], c1[0] - c1[-1]),
        min(c0[1] - c0[-1], c1[1] - c1[-1]),
        min(c0[2] - c0[-1], c1[2] - c1[-1])
    ])

    c_max = np.array([
        max(c0[0] + c0[-1], c1[0] + c1[-1]),
        max(c0[1] + c0[-1], c1[1] + c1[-1]),
        max(c0[2] + c0[-1], c1[2] + c1[-1])
    ])

    for i in range(3):
        if c_min[i] > vmax[i] or c_max[i] < vmin[i]:
            return False
    return True


@numba.njit(cache=True)
def _fiber_segment_aabb_in_sphere(c0, c1, r, center):
    c_min = np.array([
        min(c0[0] - c0[-1], c1[0] - c1[-1]),
        min(c0[1] - c0[-1], c1[1] - c1[-1]),
        min(c0[2] - c0[-1], c1[2] - c1[-1])
    ])

    c_max = np.array([
        max(c0[0] + c0[-1], c1[0] + c1[-1]),
        max(c0[1] + c0[-1], c1[1] + c1[-1]),
        max(c0[2] + c0[-1], c1[2] + c1[-1])
    ])

    dmin = 0
    for i in range(3):
        if center[i] < c_min[i]:
            dmin += (center[i] - c_min[i])**2
        elif center[i] > c_max[i]:
            dmin += (center[i] - c_max[i])**2
    return dmin <= r**2


"""
--------------------------------------------------------------------------------
-----------------------------------FIBERBUNDLE----------------------------------
--------------------------------------------------------------------------------
"""


def _convert_to_fiber_bundle(data, layers, dtype):
    """ converts data into FiberBundle"""

    if data is None:
        return [], None

    if isinstance(data, Fiber):
        return [data], Layers(layers)

    if isinstance(data, FiberBundle):
        return data[:], Layers(layers)

    if not isinstance(data, (list, tuple)):
        raise TypeError(f'data is not a list: {type(data)}')

    fiber_bundle = []
    for fiber in data:
        fiber_bundle.append(Fiber(fiber))

    return fiber_bundle, Layers(layers)


class FiberBundle:
    __is_frozen = False

    def __setattr__(self, key, value):
        if self.__is_frozen and not hasattr(self, key):
            raise TypeError('%r is a frozen class' % self)
        object.__setattr__(self, key, value)

    def __freeze(self):
        self.__is_frozen = True

    def __init__(self, data=None, layers=None, dtype=None):
        self._data, self._layers = _convert_to_fiber_bundle(data, layers, dtype)
        self.__freeze()

    def __getitem__(self, item):
        return self._data[item]

    def __setitem__(self, item, value):
        self._data[item] = Fiber(value)

    def __delitem__(self, item):
        del self._data[item]

    def __str__(self):
        return self._data.__str__()

    def __repr__(self):
        return self._data.__repr__()

    def copy(self):
        """ deep copy of class """
        return copy.deepcopy(self)

    def __iter__(self):
        return iter(self._data)

    def __next__(self):
        return next(self._data)

    def __len__(self):
        return len(self._data)

    @property
    def dtype(self):
        """ dtype of containing Fibers """
        if len(self) > 1:
            return self._data[0].dtype
        else:
            return None

    def as_array(self):
        """
        Returns copy data as list(np.array)

        Returns
        -------
        res : [(n,4)-array]
            fiber bundle as list(np.array)
        """
        return [f.as_array() for f in self]

    @property
    def layers(self):
        """
        Returns layer properties of fiber_bundle

        Returns
        -------
        res : Layers
            Layers class containing [Layer].
        """
        return self._layers

    @layers.setter
    def layers(self, value):
        self._layers = Layers(value)

    def append(self, fiber):
        """ appends Fiber to FiberBundle """
        self._data.append(Fiber(fiber))

    def extend(self, fibers):
        """ extends Fiber to FiberBundle """
        for fiber in fibers:
            self._data.append(Fiber(fiber))

    def cast(self, dtype):
        """
        Cast objects into new type

        Parameters
        ----------
        dtype : type

        Returns
        -------
        res : fiber_bundle
            casted fiber_bundle
        """

        fiber_bundle = FiberBundle()
        for fiber in self:
            fiber_bundle.append(fiber.cast(dtype))
        return fiber_bundle

    def scale(self, scale, mode='all'):
        """
        Rescales fiber_bundle

        Parameters
        ----------
        scale : float
            scale factor
        mode : str, optional
            'all', 'points' or 'radii' will be scaled

        Returns
        -------
        res : FiberBundle
            scaled fiber_bundle
        """

        fiber_bundle = FiberBundle()
        for fiber in self:
            fiber_bundle.append(fiber.scale(scale, mode))
        return fiber_bundle

    def rotate(self, rot, offset=None):
        """
        Rotates fiber_bundle around offset

        Parameters
        ----------
        rot : (3,3)-array_like
            rotation matrix
        offset : 3d-array-array_like, optional
            offset for rotation center

        Returns
        -------
        res : FiberBundle
            rotated fiber_bundle
        """

        rot = np.array(rot, copy=False)
        if offset is not None:
            offset = np.array(offset, copy=False)

        fiber_bundle = FiberBundle()
        for fiber in self:
            fiber_bundle.append(fiber.rotate(rot, offset))
        return fiber_bundle

    def translate(self, offset):
        """
        Translates fiber_bundle

        Parameters
        ----------
        offset : 3d-array-array_like
            offset to translate

        Returns
        -------
        res : FiberBundle
            translated fiber_bundle
        """

        offset = np.array(offset, copy=False)

        fiber_bundle = FiberBundle()
        for fiber in self:
            fiber_bundle.append(fiber.translate(offset))
        return fiber_bundle

    def apply(self, fun):
        """
        Applies function to fibers

        Parameters
        ----------
        fun : function

        Returns
        -------
        res : FiberBundle
            fun(fiber_bundle)
        """

        fiber_bundle = FiberBundle()
        for fiber in self:
            fiber_bundle.append(fiber.apply(fun))
        return fiber_bundle

    def apply_to_points(self, fun):
        """
        Applies function to fibers positions

        Parameters
        ----------
        fun : function

        Returns
        -------
        res : FiberBundle
            fun(fiber_bundle.points)
        """

        fiber_bundle = FiberBundle()
        for fiber in self:
            fiber_bundle.append(fiber.apply_to_points(fun))
        return fiber_bundle

    def apply_to_radii(self, fun):
        """
        Applies function to fibers radii

        Parameters
        ----------
        fun : function

        Returns
        -------
        res : FiberBundle
            fun(fiber_bundle.radii)
        """

        fiber_bundle = FiberBundle()
        for fiber in self:
            fiber_bundle.append(fiber.apply_to_radii(fun))
        return fiber_bundle

    def cut(self, voi):
        """
        Cut fiber into voi. The cutting process can create multiple fibers.
        It checks every fiber_segment_aabb if it overlapps with the voi.

        Parameters
        ----------
        voi : [xmin, ymin, zmin],[xmax,ymax,zmax]
            Volume of interest of which fibers to include. E.g. same as in
            Simulation

        Returns
        -------
        res : FiberBundle
            cutted fiber_bundle
        """

        fiber_bundle = FiberBundle()
        for fiber in self:
            fiber_bundle.extend(fiber.cut(voi))
        return fiber_bundle

    def cut_sphere(self, radius, center=(0, 0, 0)):
        """
        Cut fiber into sphere. The cutting process can create multiple fibers.
        It checks every fiber_segment_aabb if it overlapps with the sphere.

        Parameters
        ----------
        radius : float
            radius of cutting sphere
        center : 3d-array
            center of cutting sphere

        Returns
        -------
        res : FiberBundle
            cutted fiber_bundle
        """

        fiber_bundle = FiberBundle()
        for fiber in self:
            fiber_bundle.extend(fiber.cut_sphere(radius, center))
        return fiber_bundle


"""
--------------------------------------------------------------------------------
----------------------------------FIBERBUNDLES----------------------------------
--------------------------------------------------------------------------------
"""


def _convert_to_fiber_bundles(data, layers, dtype):
    """ converts data into FiberBundle"""

    if data is None:
        return []

    if isinstance(data, Fiber):
        return [FiberBundle(data, layers)]

    if isinstance(data, FiberBundle):
        return [data]

    if isinstance(data, FiberBundles):
        return data[:]

    if not isinstance(data, (list, tuple)):
        raise TypeError('data is not a list')

    if layers is not None:
        if len(data) != len(layers):
            raise TypeError('[FiberBundle] and [Layers] differ in length')
    else:
        layers = [None] * len(data)

    fiber_bundles = []
    for fiber_bundle, lys in zip(data, layers):
        fiber_bundles.append(FiberBundle(fiber_bundle, lys))

    return fiber_bundles


class FiberBundles():
    __is_frozen = False

    def __setattr__(self, key, value):
        if self.__is_frozen and not hasattr(self, key):
            raise TypeError('%r is a frozen class' % self)
        object.__setattr__(self, key, value)

    def __freeze(self):
        self.__is_frozen = True

    def __init__(self, data=None, layers=None, dtype=None):
        self._data = _convert_to_fiber_bundles(data, layers, dtype)
        self.__freeze()

    def __getitem__(self, item):
        return self._data[item]

    def __setitem__(self, item, value):
        self._data[item] = FiberBundle(value)

    def __delitem__(self, item):
        del self._data[item]

    def __str__(self):
        return self._data.__str__()

    def __repr__(self):
        return self._data.__repr__()

    def copy(self):
        """ deep copy of class """
        return copy.deepcopy(self)

    def __iter__(self):
        return iter(self._data)

    def __next__(self):
        return next(self._data)

    def __len__(self):
        return len(self._data)

    @property
    def dtype(self):
        if len(self) > 1:
            return self[0].dtype
        else:
            return None

    def as_array(self):
        """
        Returns copy data as list(list(np.array))

        Returns
        -------
        res : [[(n,4)-array]]
            fiber bundle as list(list(np.array))
        """
        return [fb.as_array() for fb in self]

    @property
    def layers(self):
        """
        Returns layer properties of fiber_bundles

        Returns
        -------
        res : [Layers]
            [Layers] class containing [Layer].
            The element position corresponds to FiberBundle index
        """
        return [fb.layers for fb in self]

    @layers.setter
    def layers(self, value):
        if len(value) != len(self):
            raise ValueError('Wrong number of [layers]')
        for fb, lys in zip(self, value):
            fb.layers = lys

    def append(self, fiber_bundle):
        """ Appends FiberBundle """
        self._data.append(FiberBundle(fiber_bundle))

    def extend(self, fiber_bundles):
        """ Extends FiberBundle """
        for fiber_bundle in fiber_bundles:
            self._data.append(FiberBundle(fiber_bundle))

    def cast(self, dtype):
        """
        Cast objects into new type

        Parameters
        ----------
        dtype : type

        Returns
        -------
        res : FiberBundles
            fiber_bundles
        """

        fiber_bundles = FiberBundles()
        for fb in self:
            fiber_bundles.append(fb.cast(dtype))
        return fiber_bundles

    def scale(self, scale, mode='all'):
        """
        Rescales fiber_bundles

        Parameters
        ----------
        scale : float
            scale factor
        mode : str, optional
            'all', 'points' or 'radii' will be scaled

        Returns
        -------
        res : FiberBundles
            scaled fiber_bundles
        """

        fiber_bundles = FiberBundles()
        for fb in self:
            fiber_bundles.append(fb.scale(scale, mode))
        return fiber_bundles

    def rotate(self, rot, offset=None):
        """
        Rotates fiber_bundles around offset

        Parameters
        ----------
        fiber_bundles : [[(,4)-array, ...]]
            list of fibers
        rot : (3,3)-array_like
            rotation matrix
        offset : 3d-array-array_like, optional
            offset for rotation center

        Returns
        -------
        res : FiberBundles
            rotated fiber_bundles
        """

        rot = np.array(rot, copy=False)
        if offset is not None:
            offset = np.array(offset, copy=False)

        fiber_bundles = FiberBundles()
        for fb in self:
            fiber_bundles.append(fb.rotate(rot, offset))
        return fiber_bundles

    def translate(self, offset):
        """
        Translates fiber_bundles

        Parameters
        ----------
        offset : 3d-array-array_like
            offset to translate

        Returns
        -------
        res : FiberBundles
            translated fiber_bundles
        """

        offset = np.array(offset, copy=False)

        fiber_bundles = FiberBundles()
        for fb in self:
            fiber_bundles.append(fb.translate(offset))
        return fiber_bundles

    def apply(self, fun):
        """
        Applies function to fibers

        Parameters
        ----------
        fun : function

        Returns
        -------
        res : FiberBundles
            fun(fiber_bundles)
        """

        fiber_bundles = FiberBundles()
        for fb in self:
            fiber_bundles.append(fb.apply(fun))
        return fiber_bundles

    def apply_to_points(self, fun):
        """
        Applies function to fibers positions

        Parameters
        ----------
        fun : function

        Returns
        -------
        res : FiberBundles
            fun(fiber_bundles[...].points)
        """

        fiber_bundles = FiberBundles()
        for fb in self:
            fiber_bundles.append(fb.apply_to_points(fun))
        return fiber_bundles

    def apply_to_radii(self, fun):
        """
        Applies function to fibers radii

        Parameters
        ----------
        fun : function

        Returns
        -------
        res : FiberBundles
            fun(fiber_bundles[...].radii)
        """

        fiber_bundles = FiberBundles()
        for fb in self:
            fiber_bundles.append(fb.apply_to_radii(fun))
        return fiber_bundles

    def cut(self, voi):
        """
        Cut fiber into voi. The cutting process can create multiple fibers.
        It checks every fiber_segment_aabb if it overlapps with the voi.

        Parameters
        ----------
        voi : [xmin, ymin, zmin],[xmax,ymax,zmax]
            Volume of interest of which fibers to include. E.g. same as in
            Simulation

        Returns
        -------
        res : FiberBundles
            cutted fiber_bundles
        """

        fiber_bundles = FiberBundles()
        for fb in self:
            fiber_bundles.append(fb.cut(voi))

        return fiber_bundles

    def cut_sphere(self, radius, center=(0, 0, 0)):
        """
        Cut fiber into voi. The cutting process can create multiple fibers.
        It checks every fiber_segment_aabb if it overlapps with the voi.

        Parameters
        ----------
        radius : float
            radius of cutting sphere
        center : 3d-array
            center of cutting sphere

        Returns
        -------
        res : FiberBundles
            cutted fiber_bundles
        """

        center = np.array(center, copy=False)

        fiber_bundles = FiberBundles()
        for fb in self:
            fiber_bundles.append(fb.cut_sphere(radius, center))

        return fiber_bundles
