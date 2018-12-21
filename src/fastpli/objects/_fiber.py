from ._fiber_cpp import _FiberCPP
import numpy as np


class Fiber(_FiberCPP):

    def __init__(self, points, radii):
        points = np.asarray(points, dtype=np.float32)
        radii = np.asarray(radii, dtype=np.float32)
        super().__init__(points, radii)
        # _FiberCPP.__init__(self, points, radii)

    @property
    def data(self):
        return (self.points, self.radii)

    def rotate_fiber(self, mat):
        mat = mat.flatten()
        super().rotate(mat)
        # _FiberCPP.rotate(self, mat)

    def rotate_around_point(self, mat, p):
        mat = mat.flatten()
        super().rotate_around_point(mat, p)
        # _FiberCPP.rotate_around_point(self, mat, p)
