import numpy as np
from ...tools import rotation
from . import fill


def cylinder(p0, p1, r0, r1, alpha, beta, mode, spacing, steps):
    p0 = np.array(p0)
    p1 = np.array(p1)

    dp = p1 - p0
    rot = rotation.a_on_b(np.array((0, 0, 1)), dp)

    if mode == 'parallel':
        points = fill.circle(r1, spacing)
        points = points[points[:, 0]**2 + points[:, 1]**2 >= r0**2, :]

        points = np.dot(rot, points.T).T

        data = []
        for p in points:
            data.append(np.array((p + p0, np.array(p + p1))))

    elif mode == 'radial':
        pass
    else:
        raise ValueError('mode hast to be "parallel" or "radial"')

    return data
