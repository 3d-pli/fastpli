from .__fiber import _Fiber
import numpy as np


class Fiber(_Fiber):

    def __init__(self, points, radii):
        super().__init__(np.asarray(points, dtype=np.float32),
                         np.asarray(radii, dtype=np.float32))

    @property
    def data(self):
        return (self.points, self.radii)

    def rotate(self, mat):
        super().rotate(mat.flatten())

    def rotate_around_point(self, mat, p):
        super().rotate_around_point(mat.flatten(), p)


def add_radius_to_fiber_bundle(fiber_bundle, radius):

    fiber_bundle_new = []
    for i in range(len(fiber_bundle)):
        fiber = np.empty((fiber_bundle[i].shape[0], 4))
        fiber[:, :-1] = fiber_bundle[i]
        fiber[:, -1] = radius
        fiber_bundle_new.append(fiber.copy())
    return fiber_bundle_new
