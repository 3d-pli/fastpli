from .__cell import __Cell
import numpy as np

class Cell(__Cell):

    def __init__(self, points, radii):
        points = np.asarray(points, dtype=np.float32)
        radii = np.asarray(radii, dtype=np.float32)
        super().__init__(points, radii)

    @property
    def data(self):
        return (self.points, self.radii)

    def rotate(self, mat):
        mat = mat.flatten()
        super().rotate(mat)

    def rotate_around_point(self, mat, p):
        mat = mat.flatten()
        super().rotate_around_point(mat, p)
