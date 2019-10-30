import numpy as np
import copy

from . import fiber


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
        for i, f in enumerate(fb):
            fiber_bundles[j][i] = fiber.Rotate(fiber_bundles[j][i], rot, offset)
    return fiber_bundles


def Translate(fiber_bundles, offset):
    fiber_bundles = copy.deepcopy(fiber_bundles)
    offset = np.array(offset, copy=False)
    for i, fb in enumerate(fiber_bundle):
        for i, f in enumerate(fb):
            fiber_bundles[j][i] = fiber.Translate(fiber_bundles[j][i], offset)
    return fiber_bundles
