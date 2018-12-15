from ._fiber_cpp import _FiberCPP
import numpy as np


class Fiber(_FiberCPP):

    def __init__(self, points, radii):
        points = np.asarray(points, dtype=np.float32)
        radii = np.asarray(radii, dtype=np.float32)
        _FiberCPP.__init__(self, points, radii)

    @property
    def data(self):
        return (self.points, self.radii)

    def rotate_fiber(self, mat):
        mat = mat.flatten()
        _FiberCPP.rotate(self, mat)

    def rotate_around_point(self, mat, p):
        mat = mat.flatten()
        _FiberCPP.rotate_around_point(self, mat, p)

    # def translate(self, offset):
    #     mat = mat.flatten()
    #     _FiberCPP.translate(self, offset)

    # def scale_points(self, f):
    #     mat = mat.flatten()
    #     _FiberCPP.scale_point(self, f)

    # def scale_radii(self, f):
    #     mat = mat.flatten()
    #     _FiberCPP.scale_radii(self, f)

    # def scale(self, f):
    #     mat = mat.flatten()
    #     _FiberCPP.scale(self, f)
