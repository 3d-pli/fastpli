import numpy as np
from numba import njit
from ...tools import rotation


@njit(cache=True)
def _rot_a_on_b(a, b):
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    if np.all(a == b):
        return np.identity(3)

    v = np.cross(a, b)
    s = np.linalg.norm(v)
    c = np.dot(a, b)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    R = np.identity(3, np.float64) + vx + np.dot(vx, vx) * (1 - c) / s**2
    return R


@njit(cache=True)
def _rot_z(phi):
    return np.array(((np.cos(phi), -np.sin(phi), 0),
                     (np.sin(phi), np.cos(phi), 0), (0, 0, 1)), np.float64)


def bundle(traj, seeds, radii, scale=1):
    """
    Generates a fiber bundle along a trajectory. The fibers will be generated coresponding to their seed points. They can be scaled along the trajectory.

    Parameters
    ----------
    traj : (n,3)-array_like
        fiber_bundle trajectory
    seeds : (m,2)-array_like
        fiber seed points on a 2d plane
    radii : float or (m,)-array_like
        fiber seeds constant or individual radii
    scale : (n,)-array_like, optional
        scale seeds along traj     
    
    Returns
    -------
    res : list((nx4)-array)
        list of fibers with (x,y,z,r)-coordinates
    """

    traj = np.array(traj, copy=False, dtype=float)
    seeds = np.array(seeds, dtype=float, ndmin=2)
    radii = np.array(radii, copy=False, dtype=float, ndmin=1)
    scale = np.array(scale, copy=False, dtype=float, ndmin=1)

    if traj.ndim != 2 or traj.shape[1] != 3:
        raise TypeError('traj : (nx3)-array')
    if traj.shape[0] < 2:
        raise ValueError('traj size is to small')

    if seeds.ndim != 2 or seeds.shape[1] != 2:
        raise ValueError("seeds : (mx2)-array")
    seeds = np.insert(seeds, 2, 0, axis=1)

    if radii.size == 1:
        radii = np.repeat(radii, seeds.shape[0])
    if radii.ndim != 1:
        raise TypeError('radii : (m)-array')
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
        R = _rot_a_on_b(tangent_old, tangent_new)

        for j in range(0, seeds.shape[0]):
            seeds[j, :] = np.dot(R, seeds[j, :])
            fiber_bundle[j][i, :] = np.append(
                seeds[j, :] * scale[i] + traj[i, :3], radii[j])

        tangent_old = tangent_new.copy()

    return fiber_bundle


@njit(cache=True)
def _cylinder_parallel(p, q, seeds, r_in, r_out, alpha, beta):
    dp = q - p
    rot = _rot_a_on_b(np.array((0, 0, 1)), dp)

    # crop seeds
    seeds = seeds[seeds[:, 0]**2 + seeds[:, 1]**2 >= r_in**2, :]
    seeds = seeds[seeds[:, 0]**2 + seeds[:, 1]**2 <= r_out**2, :]

    phi = np.arctan2(seeds[:, 1], seeds[:, 0])

    # project angles -> [0, 2*np.pi)
    phi = phi % (2.0 * np.pi)
    phi[phi < 0] += 2.0 * np.pi

    # crop seeds
    seeds = seeds[(phi > alpha) & (phi < beta), :]
    seeds = np.dot(rot, seeds.T).T

    fiber_bundle = []
    for s in seeds:
        fiber_bundle.append(np.array([s + p, s + q]))

    return fiber_bundle


@njit(cache=True)
def _cylinder_circular(p, q, seeds, r_in, r_out, alpha, beta, steps):
    # create first z-zylinder which is afterwards rotated
    dp = q - p
    b = np.linalg.norm(dp)

    # crop seeds
    seeds = seeds[seeds[:, 0] >= r_in]
    seeds = seeds[seeds[:, 0] <= r_out]
    seeds = seeds[seeds[:, 1] >= 0]
    seeds = seeds[seeds[:, 1] <= b]

    # rotate plane into first position (x-z)
    r = (r_in + r_out) / 2.0
    rot = _rot_a_on_b(np.array((0, 0, 1.0)), np.array((0, 1.0, 0)))
    seeds = np.ascontiguousarray(np.dot(rot, seeds.T).T)

    # keep rotating plane along cylinder
    rot = _rot_z((beta - alpha) / (steps - 1))
    sub_data = np.empty((steps, 3))
    fiber_bundle = []
    for j in range(seeds.shape[0]):
        for i in range(steps):
            sub_data[i, :] = seeds[j, :]
            seeds[j, :] = np.dot(rot, seeds[j, :])
        fiber_bundle.append(sub_data.copy())

    # rotate cylinder into final position
    rot = _rot_z(alpha)
    for i in range(len(fiber_bundle)):
        fiber_bundle[i] = np.ascontiguousarray(np.dot(rot, fiber_bundle[i].T).T)

    rot = _rot_a_on_b(np.array((0, 0, 1.0)), dp)
    for i in range(len(fiber_bundle)):
        fiber_bundle[i] = np.ascontiguousarray(np.dot(rot, fiber_bundle[i].T).T)

    for i in range(len(fiber_bundle)):
        fiber_bundle[i] = fiber_bundle[i] + p.T

    return fiber_bundle


@njit(cache=True)
def _cylinder_radial(p, q, seeds, r_in, r_out, alpha, beta):
    dp = q - p
    a = r_in * (beta - alpha)
    b = np.linalg.norm(dp)

    seeds = seeds[seeds[:, 0] >= 0]
    seeds = seeds[seeds[:, 0] <= a]
    seeds = seeds[seeds[:, 1] >= 0]
    seeds = seeds[seeds[:, 1] <= b]

    fiber_bundle = []

    # wrap seeds onto inner cylinder wall and extend to outer wall
    for p in seeds:
        x0 = r_in * np.cos(p[0] / a * (beta - alpha))
        y0 = r_in * np.sin(p[0] / a * (beta - alpha))
        x1 = r_out * np.cos(p[0] / a * (beta - alpha))
        y1 = r_out * np.sin(p[0] / a * (beta - alpha))

        fiber_bundle.append(np.array([[x0, y0, p[1]], [x1, y1, p[1]]]))

    # rotate cylinder into final position
    rot = _rot_z(alpha)
    for i in range(len(fiber_bundle)):
        fiber_bundle[i] = np.dot(rot, fiber_bundle[i].T).T

    rot = _rot_a_on_b(np.array((0, 0, 1.0)), dp)
    for i in range(len(fiber_bundle)):
        fiber_bundle[i] = np.dot(rot, fiber_bundle[i].T).T

    for i in range(len(fiber_bundle)):
        fiber_bundle[i] = fiber_bundle[i] + (p.T + q.T) * 0.5

    return fiber_bundle


@njit(cache=True)
def add_radii(fiber_bundle, radii):
    for i, (f, r) in enumerate(zip(fiber_bundle, radii)):
        fiber_bundle[i] = np.concatenate((f, np.ones((f.shape[0], 1)) * r),
                                         axis=1)
    return fiber_bundle


def cylinder(p,
             q,
             r_in,
             r_out,
             seeds,
             radii=1,
             alpha=0,
             beta=2 * np.pi,
             mode='p',
             steps=360):
    """
    Generated a fiber bundle inside a cylinder

    Parameters
    ----------
    p,q : (3,)-array_like
        (x,y,z)-points of begin and end of cylinder
    r_in,r_out : 
        inner and outer radius of cylinder
    seeds : (m,2)-array_like
        fiber seeds on a 2d plane
    radii : float or (n,)-array_like
        fiber seeds radii
    alpha, beta : float
        fibers are between alpha end beta inside cylinder
    mode : char or string
        \"p\", \"parallel\", \"c\", \"circular\", \"r\", \"radial\"
    steps : int
        steps along fibers in circular mode   
    
    Returns
    -------
    res : list((nx4)-array)
        list of fibers with (x,y,z,r)-coordinates
    """

    p = np.array(p, dtype=float)
    q = np.array(q, dtype=float)
    r_in = float(r_in)
    r_out = float(r_out)
    seeds = np.array(seeds, dtype=float, ndmin=2)
    radii = np.array(radii, dtype=float, ndmin=1)
    alpha = float(alpha)
    beta = float(beta)
    steps = int(steps)

    if alpha == beta:
        return np.array([])

    if seeds.ndim != 2 or seeds.shape[1] != 2:
        raise ValueError("seeds : (nx2)-array")
    seeds = np.insert(seeds, 2, 0, axis=1)

    if radii.ndim != 1:
        raise TypeError('radii : (n)-array')
    if radii.size == 1:
        radii = np.repeat(radii, seeds.shape[0])
    if radii.size != seeds.shape[0]:
        raise ValueError('radii must have the same length as seeds')

    dp = q - p
    # fiber_bundle = []
    if mode == 'parallel' or mode == 'p':
        fiber_bundle = _cylinder_parallel(p, q, seeds, r_in, r_out, alpha, beta)

    elif mode == 'circular' or mode == 'c':
        fiber_bundle = _cylinder_circular(p, q, seeds, r_in, r_out, alpha, beta,
                                          steps)

    elif mode == 'radial' or mode == 'r':
        fiber_bundle = _cylinder_radial(p, q, seeds, r_in, r_out, alpha, beta)

    else:
        raise ValueError('mode has to be "parallel" or "radial"')

    add_radii(fiber_bundle, radii)

    return fiber_bundle


@njit(cache=True)
def _ray_box_intersection(p, dir, b_min, b_max):

    tmin = np.divide(b_min - p, dir)
    tmax = np.divide(b_max - p, dir)

    tmin, tmax = np.minimum(tmin, tmax), np.maximum(tmin, tmax)

    return np.max(tmin), np.min(tmax)


@njit(cache=True)
def _ray_box_intersection_pp(p, q, b_min, b_max):
    return _ray_box_intersection(p, q - p, b_min, b_max)


def cuboid(p, q, phi, theta, seeds, radii):
    """
    Generated a fiber bundle inside a cuboid along a axis

    Parameters
    ----------
    p,q : (3,)-array_like
        (x,y,z)-points of cuboid corners
    phi, theta : float
        spherical angles of axis
    seeds : (m,2)-array_like
        fiber seeds on a 2d plane
    radii : float or (n,)-array_like
        fiber seeds radii
    
    Returns
    -------
    res : list((nx4)-array)
        list of fibers with (x,y,z,r)-coordinates
    """

    p = np.array(p, copy=False)
    q = np.array(q, copy=False)
    seeds = np.array(seeds, float)
    phi = float(phi)
    theta = float(theta)
    radii = np.array(radii, ndmin=1)

    if p.ndim != 1 or q.ndim != 1:
        raise ValueError("p.ndim, q.ndim = 1")
    if p.size != 3 or q.size != 3:
        raise ValueError("p.size, q.size != 3")

    if seeds.ndim != 2 or seeds.shape[1] != 2:
        raise ValueError("seeds : (nx2)-ndarray")
    seeds = np.insert(seeds, 2, 0, axis=1)

    if radii.shape[0] == 1:
        radii = np.ones(seeds.shape[0], radii.dtype) * radii

    if seeds.shape[0] != radii.shape[0]:
        raise ValueError("seeds.shape[1] != radii.shape[0]")

    # rotate fibers
    dir = np.array([
        np.cos(phi) * np.sin(theta),
        np.sin(phi) * np.sin(theta),
        np.cos(theta)
    ])
    rot = _rot_a_on_b(np.array([0, 0, 1.0]), dir)
    seeds = np.dot(rot, seeds.T).T + 0.5 * (p + q)

    fiber_bundle = []

    for s, r in zip(seeds, radii):
        # find ray box intersection
        s -= 0.5 * dir
        t_min, t_max = _ray_box_intersection_pp(s, s + dir, p, q)

        if t_min >= t_max:  # outside of volume
            continue

        p_min = s + t_min * dir
        p_max = s + t_max * dir

        fiber_bundle.append(
            np.array([[p_min[0], p_min[1], p_min[2], r],
                      [p_max[0], p_max[1], p_max[2], r]]))

    return fiber_bundle
