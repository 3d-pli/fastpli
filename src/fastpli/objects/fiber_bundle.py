import numpy as np
import copy

from . import fiber


def Rescale(fiber_bundle, scale, mod='all'):
    fiber_bundle = copy.deepcopy(fiber_bundle)
    for i, f in enumerate(fiber_bundle):
        fiber_bundle[i] = fiber.Rescale(fiber_bundle[i], scale, mod)
    return fiber_bundle


def Rotate(fiber_bundle, rot, offset=None):
    fiber_bundle = copy.deepcopy(fiber_bundle)
    rot = np.array(rot, copy=False)
    if offset is not None:
        offset = np.array(offset, copy=False)
    for i, f in enumerate(fiber_bundle):
        fiber_bundle[i] = fiber.Rotate(fiber_bundle[i], rot, offset)
    return fiber_bundle


def Translate(fiber_bundle, offset):
    fiber_bundle = copy.deepcopy(fiber_bundle)
    offset = np.array(offset, copy=False)
    for i, f in enumerate(fiber_bundle):
        fiber_bundle[i] = fiber.Translate(fiber_bundle[i], offset)
    return fiber_bundle
