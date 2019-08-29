import numpy as np
from ...tools.rotation import a_on_b


def _hexgrid(a, b, spacing):
    # dx = np.sin(np.deg2rad(30)) * spacing
    # dy = np.cos(np.deg2rad(30)) * spacing

    dx = 0.5 * spacing
    dy = 0.8660254037844386 * spacing

    points_0 = np.mgrid[0:b + spacing / 2:2 * dy, 0:a +
                        spacing / 2:spacing].reshape(2, -1)[1::-1]
    points_1 = np.mgrid[dy:b + spacing / 2:2 * dy, dx:spacing / 2 +
                        a:spacing].reshape(2, -1)[1::-1]

    return np.concatenate((points_0, points_1), axis=1).T


def rectangle(a, b, spacing, mode='center'):

    if mode == 'center':
        offset = np.array((a / 2.0, b / 2.0))
    elif mode == 'origin':
        offset = np.array((0, 0))
    else:
        ValueError('mode has to be "center" or "origin"')

    return _hexgrid(a, b, spacing) - offset


def circle(radius, spacing):

    points = _hexgrid(2 * radius, 2 * radius, spacing) - np.array(
        (radius, radius))
    dr2 = points[:, 0]**2 + points[:, 1]**2
    points = points[dr2 <= radius**2, :]

    return points


def bundle(traj, seeds, fiber_radius):
    '''
    traj: fiber_bundle trajectory: (nx3 or nx4)-array
    object: fiber seed on a 2d plane: (nx2 or nx3)-array
    fiber_radius: (int, float, (nx1)-array)
    '''

    traj = np.array(traj, dtype=np.float64, copy=False)
    seeds = np.atleast_2d(np.array(seeds, dtype=np.float64, copy=True))

    if len(traj.shape) != 2 or not (traj.shape[1] != 3 or traj.shape[1] != 4):
        raise TypeError('traj format: (nx3)-array or (nx4)-array')

    if traj.shape[1] == 3:
        traj = np.append(traj, np.ones((traj.shape[0], 1)), axis=1)

    if traj.shape[0] < 2:
        raise TypeError('traj size is to small')

    if len(seeds.shape) != 2 or not (seeds.shape[1] != 2 or
                                     seeds.shape[1] != 3):
        raise TypeError('seeds format: (nx2)-array or (nx3)-array')

    if seeds.shape[1] == 2:
        seeds = np.append(seeds, np.zeros((seeds.shape[0], 1)), axis=1)

    if not isinstance(fiber_radius, (int, float, np.ndarray)):
        raise TypeError('fiber_radius format: (int, float, np.ndarray)')

    if isinstance(fiber_radius, (int, float)):
        fiber_radius = np.ones(seeds.shape[0]) * fiber_radius

    if len(fiber_radius.shape) != 1 or fiber_radius.size != seeds.shape[0]:
        raise TypeError('fiber_radius must have the same length as seeds')

    fiber_radius = np.array(fiber_radius, dtype=np.float64, copy=False)
    fiber_bundle = [np.empty([traj.shape[0], 4]) for i in range(seeds.shape[0])]
    tangent_old = np.array([0, 0, 1], float)

    for i in range(0, traj.shape[0]):
        if i == 0:
            tangent_new = 0.5 * (traj[i + 1, :3] - traj[i, :3])
        elif i == traj.shape[0] - 1:
            tangent_new = 0.5 * (traj[i, :3] - traj[i - 1, :3])
        else:
            tangent_new = 0.5 * (traj[i + 1, :3] - traj[i - 1, :3])

        if np.all(tangent_new == 0):
            raise TypeError('same point:', i)

        tangent_new = tangent_new / np.linalg.norm(tangent_new)
        R = a_on_b(tangent_old, tangent_new)

        for j in range(0, seeds.shape[0]):
            seeds[j, :] = np.dot(R, seeds[j, :])
            fiber_bundle[j][i, :] = np.append(
                seeds[j, :] * traj[i, -1] + traj[i, :3], fiber_radius[j])

        tangent_old = tangent_new.copy()

    return fiber_bundle
