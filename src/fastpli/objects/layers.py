"""
layer classes
"""

import warnings
import copy


def _check_layer(ly):
    if len(ly) != 4:
        raise TypeError('layer != (float, float, float, char)')
    if ly[1] < 0 and ly[-1] == 'r':
        warnings.warn('birefringence negative and radial')
    if ly[1] > 0 and ly[-1] == 'p':
        warnings.warn('birefringence positive and parallel')
    if ly[1] != 0 and ly[-1] == 'b':
        raise ValueError('birefringence != 0 for background.')
    return (float(ly[0]), float(ly[1]), float(ly[2]), ly[3])


class Layer:
    __is_frozen = False

    def __setattr__(self, key, value):
        if self.__is_frozen and not hasattr(self, key):
            raise TypeError('%r is a frozen class' % self)
        object.__setattr__(self, key, value)

    def __freeze(self):
        self.__is_frozen = True

    def __init__(self, *args):
        if len(args) == 1:
            args = args[0]
        self._data = _check_layer(args)
        self.__freeze()

    def __getitem__(self, item):
        return self._data[item]

    def __str__(self):
        return self._data.__str__()

    def __repr__(self):
        return self._data.__repr__()

    def copy(self):
        return copy.deepcopy(self)

    def data(self):
        return self._data

    @property
    def scale(self):
        return self._data[0]

    @property
    def birefringence(self):
        return self._data[1]

    @property
    def absorption(self):
        return self._data[2]

    @property
    def model(self):
        return self._data[3]


def _check_layers(lys):
    if lys is None:
        return None

    if all(isinstance(ly, (list, tuple, Layer)) for ly in lys):
        lys = [Layer(ly) for ly in lys]

    elif isinstance(lys, (list, tuple, Layer)):
        lys = [Layer(lys)]
    else:
        raise TypeError('Expected \"Layer\" or \"[Layer]\"')

    if any(ly.scale > 1 for ly in lys):
        warnings.warn('layer scale > 1 detected')

    keys = sorted(range(len(lys)), key=lambda i: lys[i][0])
    if any(keys[i] > keys[i + 1] for i in range(len(keys) - 1)):
        warnings.warn('layers had to be sorted')
        lys = [lys[i] for i in keys]

    if any(lys[i][0] == lys[i + 1][0] for i in range(len(lys) - 1)):
        raise ValueError('detected layers of same radii')

    return lys


class Layers:
    __is_frozen = False

    def __setattr__(self, key, value):
        if self.__is_frozen and not hasattr(self, key):
            raise TypeError('%r is a frozen class' % self)
        object.__setattr__(self, key, value)

    def __freeze(self):
        self.__is_frozen = True

    def __init__(self, data):
        self._data = _check_layers(data)
        self.__freeze()

    def __getitem__(self, item):
        return self._data[item]

    def __setitem__(self, item, value):
        lys = self._data
        lys[item] = Layer(value)

        keys = sorted(range(len(lys)), key=lambda i: lys[i][0])
        if any(keys[i] > keys[i + 1] for i in range(len(keys) - 1)):
            warnings.warn('layers had to be sorted')
            lys = [lys[i] for i in keys]

        if any(lys[i][0] == lys[i + 1][0] for i in range(len(lys) - 1)):
            raise ValueError('detected layers of same radii')

        self._data = lys

    def __str__(self):
        return self._data.__str__()

    def __repr__(self):
        return self._data.__repr__()

    def copy(self):
        return copy.deepcopy(self)

    def __iter__(self):
        return self._data.__iter__()

    def __next__(self):
        return self._data.__next__()

    def __len__(self):
        return self._data.__len__()

    def data(self):
        return [ly.data() for ly in self]
