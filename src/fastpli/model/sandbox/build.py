import numpy as np
from ...tools import rotation


def bundle(traj, seeds, radii, scale=1):
    '''
    traj: fiber_bundle trajectory: (nx3)-array
    seeds: fiber seeds on a 2d plane: (mx2 or mx3)-array
    radii: fiber seeds radii: (int, float, (m)-array)
    scale: scale seeds along traj: (n)-array
    '''

    traj = np.array(traj, copy=False, dtype=float)
    seeds = np.array(seeds, dtype=float, ndmin=2)
    radii = np.array(radii, copy=False, dtype=float, ndmin=1)
    scale = np.array(scale, copy=False, dtype=float, ndmin=1)

    if traj.ndim != 2 or traj.shape[1] != 3:
        raise TypeError('traj : (nx3)-array')
    if traj.shape[0] < 2:
        raise ValueError('traj size is to small')

    if seeds.ndim != 2 or (seeds.shape[1] != 2 and seeds.shape[1] != 3):
        raise TypeError('seeds : (nx2)-array or (nx3)-array')
    if seeds.shape[1] == 2:
        seeds = np.append(seeds, np.zeros((seeds.shape[0], 1)), axis=1)

    if radii.size == 1:
        radii = np.repeat(radii, seeds.shape[0])
    if radii.ndim != 1:
        raise TypeError('radii : (n)-array')
    if radii.size != seeds.shape[0]:
        raise ValueError('radii must have the same length as seeds')

    if scale.size == 1:
        scale = np.ones(traj.shape[0]) * scale[0]
    if scale.ndim != 1:
        raise TypeError('scale : (n)-array')
    if scale.size != traj.shape[0]:
        raise ValueError('scale must have the same length as traj')

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
        R = rotation.a_on_b(tangent_old, tangent_new)

        for j in range(0, seeds.shape[0]):
            seeds[j, :] = np.dot(R, seeds[j, :])
            fiber_bundle[j][i, :] = np.append(
                seeds[j, :] * scale[i] + traj[i, :3], radii[j])

        tangent_old = tangent_new.copy()

    return fiber_bundle


def cylinder(p0,
             p1,
             r0,
             r1,
             seeds,
             radii=1,
             alpha=0,
             beta=2 * np.pi,
             mode='p',
             steps=360):
    '''
    p0,p1: (x,y,z)-points of begin and end of cylinder
    r0,r1: radius of cylinder at position 0 or 1
    radii: fiber seeds radii: (int, float, (m)-array)
    scale: scale seeds along traj: (n)-array
    alpha, beta: fibers are between alpha end beta inside cylinder
    mode: parallel, circular, radial
    steps: only for circular mode
    '''

    p0 = np.array(p0, dtype=float)
    p1 = np.array(p1, dtype=float)
    r0 = float(r0)
    r1 = float(r1)
    seeds = np.array(seeds, dtype=float, ndmin=2)
    radii = np.array(radii, dtype=float, ndmin=1)
    alpha = float(alpha)
    beta = float(beta)
    steps = int(steps)

    # project angles -> [0, 2*np.pi)

    alpha = alpha % (2.0 * np.pi)
    if alpha < 0:
        alpha += 2.0 * np.pi

    if beta != 2 * np.pi:
        beta = beta % (2.0 * np.pi)
    if beta < 0:
        beta += 2.0 * np.pi

    if seeds.ndim != 2 or seeds.shape[1] != 2:
        raise TypeError('seeds : (nx2)-array')
    if seeds.shape[1] == 2:
        seeds = np.insert(seeds, 2, 0, axis=1)

    if radii.ndim != 1:
        raise TypeError('radii : (n)-array')
    if radii.size == 1:
        radii = np.repeat(radii, seeds.shape[0])
    if radii.size != seeds.shape[0]:
        raise ValueError('radii must have the same length as seeds')

    dp = p1 - p0
    fiber_bundle = []
    if mode == 'parallel' or mode == 'p':
        rot = rotation.a_on_b(np.array((0, 0, 1)), dp)

        # crop seeds
        seeds = seeds[seeds[:, 0]**2 + seeds[:, 1]**2 >= r0**2, :]
        seeds = seeds[seeds[:, 0]**2 + seeds[:, 1]**2 <= r1**2, :]

        phi = np.arctan2(seeds[:, 1], seeds[:, 0])

        # project angles -> [0, 2*np.pi)
        phi = phi % (2.0 * np.pi)
        phi[phi < 0] += 2.0 * np.pi

        # crop seeds
        seeds = seeds[(phi > alpha) & (phi < beta), :]
        seeds = np.dot(rot, seeds.T).T

        for p, r in zip(seeds, radii):
            fiber_bundle.append(
                np.insert(np.array([p + p0, p + p1]), 3, r, axis=1))

    elif mode == 'circular' or mode == 'c':
        # create first z-zylinder which is afterwards rotated
        # a = r1 - r0
        b = np.linalg.norm(dp)

        # crop seeds
        seeds = seeds[seeds[:, 0] >= r0]
        seeds = seeds[seeds[:, 0] <= r1]
        seeds = seeds[seeds[:, 2] >= 0]
        seeds = seeds[seeds[:, 2] <= b]

        # rotate plane into first position (x-z)
        r = (r0 + r1) / 2.0
        rot = rotation.a_on_b(np.array((0, 0, 1)), np.array((0, 1, 0)))
        seeds = np.dot(rot, seeds.T).T

        # keep rotating plane along cylinder
        rot = rotation.z((beta - alpha) / (steps - 1))
        sub_data = np.empty((steps, 3))
        for p in seeds:
            for i in range(steps):
                sub_data[i, :] = p
                p = np.dot(rot, p)
            fiber_bundle.append(sub_data.copy())

        # rotate cylinder into final position
        rot = rotation.z(alpha)
        for i in range(len(fiber_bundle)):
            fiber_bundle[i] = np.dot(rot, fiber_bundle[i].T).T

        rot = rotation.a_on_b(np.array((0, 0, 1)), dp)
        for i in range(len(fiber_bundle)):
            fiber_bundle[i] = np.dot(rot, fiber_bundle[i].T).T

        for i in range(len(fiber_bundle)):
            fiber_bundle[i] = fiber_bundle[i] + p0.T

    elif mode == 'radial' or mode == 'r':
        a = r0 * (beta - alpha)
        b = np.linalg.norm(dp)

        seeds = seeds[seeds[:, 0] >= 0]
        seeds = seeds[seeds[:, 0] <= a]
        seeds = seeds[seeds[:, 1] >= 0]
        seeds = seeds[seeds[:, 1] <= b]

        # wrap seeds onto inner cylinder wall and extend to outer wall
        for p in seeds:
            x0 = r0 * np.cos(p[0] / a * (beta - alpha))
            y0 = r0 * np.sin(p[0] / a * (beta - alpha))
            x1 = r1 * np.cos(p[0] / a * (beta - alpha))
            y1 = r1 * np.sin(p[0] / a * (beta - alpha))

            fiber_bundle.append(np.array([[x0, y0, p[1]], [x1, y1, p[1]]]))

        # rotate cylinder into final position
        rot = rotation.z(alpha)
        for i in range(len(fiber_bundle)):
            fiber_bundle[i] = np.dot(rot, fiber_bundle[i].T).T

        rot = rotation.a_on_b(np.array((0, 0, 1)), dp)
        for i in range(len(fiber_bundle)):
            fiber_bundle[i] = np.dot(rot, fiber_bundle[i].T).T

        for i in range(len(fiber_bundle)):
            fiber_bundle[i] = fiber_bundle[i] + (p0.T + p1.T) * 0.5

    else:
        raise ValueError('mode has to be "parallel" or "radial"')

    return fiber_bundle


def _ray_box_intersection(p, dir, b_min, b_max):

    tmin = np.divide(b_min - p, dir)
    tmax = np.divide(b_max - p, dir)

    tmin, tmax = np.minimum(tmin, tmax), np.maximum(tmin, tmax)

    return np.max(tmin), np.min(tmax)


def _ray_box_intersection_pp(p0, p1, b_min, b_max):

    if np.any(np.maximum(p0, p1) < np.minimum(b_min, b_max)) or np.any(
            np.minimum(p0, p1) > np.maximum(b_min, b_max)):
        return 0, 0

    return _ray_box_intersection(p0, p1 - p0, b_min, b_max)


def cuboid(a, b, phi, theta, seeds, radii):

    a = np.array(a, float)
    b = np.array(b, float)
    phi = float(phi)
    theta = float(theta)
    radii = np.array(radii, ndmin=1)

    if a.ndim != 1 or b.ndim != 1:
        raise ValueError("a.ndim, b.ndim = 1")
    if a.size != b.size:
        raise ValueError("a.size != b.size")

    if seeds.ndim != 2:
        raise ValueError("")
    if seeds.shape[1] == 1 or seeds.shape > 3:
        raise ValueError("")
    if seeds.shape[1] == 2:
        seeds = np.insert(seeds, 2, 0, axis=1)

    # rotate fibers
    dir = np.array([
        np.cos(phi) * np.sin(theta),
        np.sin(phi) * np.sin(theta),
        np.cos(theta)
    ])
    rot = rotation.a_on_b(np.array([0, 0, 1]), dir)
    seeds = np.dot(rot, seeds.T).T + 0.5 * (a + b)

    fiber_bundle = []
    # dir = dir * cube_length
    for p, r in zip(seeds, radii):
        # find ray box intersection
        p -= 0.5 * dir
        t_min, t_max = _ray_box_intersection_pp(p, p + dir, a, b)
        if t_min >= t_max:  # outside of volume
            continue

        p_min = p + t_min * dir
        p_max = p + t_max * dir

        fiber_bundle.append(
            np.array([[p_min[0], p_min[1], p_min[2],
                       r][p_max[0], p_max[1], p_max[2], r]]))

    return fiber_bundle
