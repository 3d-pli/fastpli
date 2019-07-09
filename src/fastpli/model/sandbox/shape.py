import numpy as np
from ...tools import rotation
from . import fill


def cylinder(p0, p1, r0, r1, alpha, beta, mode, spacing, steps):
    # project angles -> [0, 2*np.pi)
    alpha = alpha % (2.0 * np.pi)
    if alpha < 0:
        alpha += 2.0 * np.pi
    beta = beta % (2.0 * np.pi)
    if beta < 0:
        beta += 2.0 * np.pi

    p0 = np.array(p0)
    p1 = np.array(p1)
    dp = p1 - p0
    steps = max(2, steps)

    data = []
    if mode == 'parallel' or mode == 'p':
        rot = rotation.a_on_b(np.array((0, 0, 1)), dp)

        points = fill.circle(r1, spacing)
        points = points[points[:, 0]**2 + points[:, 1]**2 >= r0**2, :]
        phi = np.arctan2(points[:, 1], points[:, 0])

        # project angles -> [0, 2*np.pi)
        phi = phi % (2.0 * np.pi)
        phi[phi < 0] = phi[phi < 0] + 2.0 * np.pi

        points = points[(phi > alpha) & (phi < beta), :]
        points = np.dot(rot, points.T).T

        for p in points:
            x = np.interp(np.arange(steps), [0, steps - 1], [p0[0], p1[0]])
            y = np.interp(np.arange(steps), [0, steps - 1], [p0[1], p1[1]])
            z = np.interp(np.arange(steps), [0, steps - 1], [p0[2], p1[2]])
            data.append(np.array([x + p[0], y + p[1], z + p[2]]).T)
    elif mode == 'circular' or mode == 'c':
        # create first z-zylinder which is afterwards rotated
        a = r1 - r0
        b = np.linalg.norm(dp)
        points = fill.rectangle(a, b, spacing, 'center')

        # rotate plane into first position
        r = (r0 + r1) / 2.0
        rot = rotation.a_on_b(np.array((0, 0, 1)), np.array((0, 1, 0)))

        # print(rot)
        points = np.dot(rot, points.T).T
        points[:, 0] += r

        # keep rotating plane along cylinder
        rot = rotation.z((beta - alpha) / (steps - 1))
        sub_data = np.empty((steps, 3))
        for p in points:
            for i in range(steps):
                sub_data[i, :] = p
                p = np.dot(rot, p)
            data.append(sub_data.copy())

        # rotate cylinder into final position
        rot = rotation.z(alpha)
        for i in range(len(data)):
            data[i] = np.dot(rot, data[i].T).T

        rot = rotation.a_on_b(np.array((0, 0, 1)), dp)
        for i in range(len(data)):
            data[i] = np.dot(rot, data[i].T).T

        for i in range(len(data)):
            data[i] = data[i] + (p0.T + p1.T) * 0.5

    elif mode == 'radial' or mode == 'r':
        a = r0 * (beta - alpha)
        b = np.linalg.norm(dp)
        points = fill.rectangle(a, b, spacing, 'center')
        points[:, 0] += a / 2

        for p in points:
            x0 = r0 * np.cos(p[0] / a * (beta - alpha))
            y0 = r0 * np.sin(p[0] / a * (beta - alpha))
            x1 = r1 * np.cos(p[0] / a * (beta - alpha))
            y1 = r1 * np.sin(p[0] / a * (beta - alpha))

            x = np.interp(np.arange(steps), [0, steps - 1], [x0, x1])
            y = np.interp(np.arange(steps), [0, steps - 1], [y0, y1])
            z = np.ones(x.size) * p[1]

            data.append(np.array([x, y, z]).T)

        # rotate cylinder into final position
        rot = rotation.z(alpha)
        for i in range(len(data)):
            data[i] = np.dot(rot, data[i].T).T

        rot = rotation.a_on_b(np.array((0, 0, 1)), dp)
        for i in range(len(data)):
            data[i] = np.dot(rot, data[i].T).T

        for i in range(len(data)):
            data[i] = data[i] + (p0.T + p1.T) * 0.5

    else:
        raise ValueError('mode has to be "parallel" or "radial"')

    return data
