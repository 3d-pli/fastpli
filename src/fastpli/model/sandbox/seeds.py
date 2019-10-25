import numpy as np


def triangle_grid(a, b, spacing):
    # dx = np.sin(np.deg2rad(30)) * spacing
    # dy = np.cos(np.deg2rad(30)) * spacing

    a = float(a)
    b = float(b)
    spacing = float(spacing)

    dx = 0.5 * spacing
    dy = 0.8660254037844386 * spacing

    points_0 = np.mgrid[0:a + spacing / 2:spacing, 0:b + spacing / 2:2 *
                        dy].reshape(2, -1).T
    points_1 = np.mgrid[dx:spacing / 2 + a:spacing, dy:b + spacing / 2:2 *
                        dy].reshape(2, -1).T

    return np.concatenate((points_0, points_1))


def crop_rectangle(a, b, seeds, radii=0):
    '''
    a,b : number = length of crop_rectangle
    a,b : (2,1)ndarray = (xmin, ymin),(xmax, ymax)
    seeds: list of 2d coordinates
    '''

    a = np.array(a, ndmin=1)
    b = np.array(b, ndmin=1)
    seeds = np.array(seeds, copy=False, ndmin=2)
    radii = np.array(radii, ndmin=1)

    if a.ndim != 1 or b.ndim != 1:
        raise ValueError("a.ndim, b.ndim != 1")
    if a.size != b.size:
        raise ValueError("a.size != b.size")

    if seeds.ndim != 2 or seeds.shape[1] != 2:
        raise TypeError('seeds : (nx2)-array')

    if a.size == 1:
        return seeds[np.logical_and(
            np.logical_and(seeds[:, 0] - radii >= 0, seeds[:, 0] + radii <= a),
            np.logical_and(seeds[:, 1] - radii >= 0, seeds[:, 1] + radii <= b))]
    elif a.size == 2:
        return seeds[np.logical_and(np.all(seeds >= a), np.all(seeds <= b))]

    else:
        raise ValueError("a.size, b.size > 2")


def crop_circle(radius, seeds, radii=0, center=[0, 0]):

    radius = float(radius)
    seeds = np.array(
        seeds,
        dtype=float,
        copy=False,
    )
    radii = np.array(radii, dtype=float, copy=False, ndmin=1)
    center = np.array(center)

    if seeds.ndim != 2 or seeds.shape[1] != 2:
        raise TypeError('seeds : (nx2)-array')

    if center.ndim != 1 or center.size != 2:
        raise ValueError("center : [x,y]")

    dr2 = (seeds[:, 0] - center[0])**2 + (seeds[:, 1] - center[1])**2

    return seeds[(dr2 + radii**2) <= radius**2, :]
