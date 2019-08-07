import numpy as np

from . import fiber


def Resize(fiber_bundle, scale, mod='all'):
    for i, f in enumerate(fiber_bundle):
        fiber_bundle[i] = fiber.Resize(fiber_bundle[i], scale, mod)
    return fiber_bundle


def Rotate(fiber_bundle, rot, offset=None):
    rot = np.array(rot, copy=False)
    if offset is not None:
        offset = np.array(offset, copy=False)
    for i, f in enumerate(fiber_bundle):
        fiber_bundle[i] = fiber.Rotate(fiber_bundle[i], rot, offset)
    return fiber_bundle


def Translate(fiber_bundle, offset):
    offset = np.array(offset, copy=False)
    for i, f in enumerate(fiber_bundle):
        fiber_bundle[i] = fiber.Translate(fiber_bundle[i], offset)
    return fiber_bundle
