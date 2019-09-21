import numpy as np

from . import fiber


def Resize(fiber_bundles, scale, mod='all'):
    for j, fb in enumerate(fiber_bundles):
        for i, f in enumerate(fb):
            fiber_bundles[j][i] = fiber.Resize(fiber_bundles[j][i], scale, mod)
    return fiber_bundles


def Rotate(fiber_bundles, rot, offset=None):
    rot = np.array(rot, copy=False)
    if offset is not None:
        offset = np.array(offset, copy=False)
    for j, fb in enumerate(fiber_bundles):
        for i, f in enumerate(fb):
            fiber_bundles[j][i] = fiber.Rotate(fiber_bundles[j][i], rot, offset)
    return fiber_bundles


def Translate(fiber_bundles, offset):
    offset = np.array(offset, copy=False)
    for i, fb in enumerate(fiber_bundle):
        for i, f in enumerate(fb):
            fiber_bundles[j][i] = fiber.Translate(fiber_bundles[j][i], offset)
    return fiber_bundles
