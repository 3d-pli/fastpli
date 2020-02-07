import numpy as np
import copy

from . import fiber


def Cast(fiber_bundles):

    if not fiber_bundles:
        return fiber_bundles

    if not isinstance(fiber_bundles, (list, tuple)):
        raise TypeError("fiber_bundles is not a list")

    for fb_i, fb in enumerate(fiber_bundles):
        if not isinstance(fb, (list, tuple)):
            raise TypeError("fiber_bundle is not a list")

        for f_i, f in enumerate(fb):
            fiber_bundles[fb_i][f_i] = np.array(f, dtype=float)

            if fiber_bundles[fb_i][f_i].ndim != 2 or fiber_bundles[fb_i][
                    f_i].shape[1] != 4:
                raise TypeError("fiber elements has to be of shape nx4")

    return fiber_bundles


def Rescale(fiber_bundles, scale, mod='all'):
    fiber_bundles = copy.deepcopy(fiber_bundles)
    for j, fb in enumerate(fiber_bundles):
        for i, f in enumerate(fb):
            fiber_bundles[j][i] = fiber.Rescale(fiber_bundles[j][i], scale, mod)
    return fiber_bundles


def Rotate(fiber_bundles, rot, offset=None):
    fiber_bundles = copy.deepcopy(fiber_bundles)
    rot = np.array(rot, copy=False)
    if offset is not None:
        offset = np.array(offset, copy=False)
    for j, fb in enumerate(fiber_bundles):
        for i, _ in enumerate(fb):
            fiber_bundles[j][i] = fiber.Rotate(fiber_bundles[j][i], rot, offset)
    return fiber_bundles


def Translate(fiber_bundles, offset):
    fiber_bundles = copy.deepcopy(fiber_bundles)
    offset = np.array(offset, copy=False)
    for i, fb in enumerate(fiber_bundles):
        for i, _ in enumerate(fb):
            fiber_bundles[j][i] = fiber.Translate(fiber_bundles[j][i], offset)
    return fiber_bundles
