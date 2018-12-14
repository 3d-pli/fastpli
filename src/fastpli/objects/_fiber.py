from ._fiber_cpp import _FiberCPP
import numpy as np


class Fiber(_FiberCPP):

    def __init__(self, points, radii):
        points = np.asarray(points, dtype=np.float32)
        radii = np.asarray(radii, dtype=np.float32)
        _FiberCPP.__init__(self, points, radii)
