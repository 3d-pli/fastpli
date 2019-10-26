import numpy as np


def triangular_grid(a, b, spacing, center=False, sort=True):
    ''' 
    2d triangular grid from [0,a)x[0,b)
    center: grid from [-a/2,a/2)x[-b/2,b/2) with point at (0,0)
    sort: lexsort x,y
    '''

    a = float(a)
    b = float(b)
    spacing = float(spacing)

    dx = spacing / 2  # np.sin(np.deg2rad(30))
    dy = spacing * np.sqrt(3) / 2  # np.cos(np.deg2rad(30))

    x0 = 0
    y0 = 0

    if center:
        x0 = -((a / 2) // (2 * dx)) * 2 * dx
        y0 = -((b / 2) // (2 * dy)) * 2 * dy
        a = a / 2
        b = b / 2

    grid_0 = np.mgrid[x0:a:spacing, y0:b:2 * dy].reshape(2, -1)
    grid_1 = np.mgrid[x0 + dx:a:spacing, y0 + dy:b:2 * dy].reshape(2, -1)
    grid = np.concatenate((grid_0, grid_1), axis=1)

    if sort:
        idx = np.lexsort((grid[0, :], grid[1, :]))
        grid = grid[:, idx]

    return np.ascontiguousarray(grid.T)


if __name__ == "__main__":
    print(triangular_grid(4.1, 6.1, 1, True))


def crop_rectangle(a, b, seeds, radii=0):
    '''
    a,b : number: (0,0)--(a,b)
    a,b : array: (xmin, ymin)--(xmax, ymax)
    seeds: list of 2d coordinates
    '''

    a = np.array(a, ndmin=1)
    b = np.array(b, ndmin=1)
    seeds = np.array(seeds, ndmin=2, copy=False)
    radii = np.array(radii, cope=False)

    if a.ndim != 1 or b.ndim != 1:
        raise ValueError("a.ndim, b.ndim != 1")
    if a.size != b.size:
        raise ValueError("a.size != b.size")
    if a.size == 1:
        a, b = np.array((0, 0)), np.array((a, b))
    elif a.size > 2:
        raise ValueError("a.size, b.size > 2")

    if seeds.ndim != 2 or seeds.shape[1] != 2:
        raise TypeError('seeds : (nx2)-array')

    return seeds[np.logical_and(np.all(seeds >= a + radius),
                                np.all(seeds <= b - radius))]


def crop_circle(radius, seeds, center=[0, 0], radii=0):

    radius = float(radius)
    seeds = np.array(seeds, ndmin=2, copy=False)
    center = np.array(center)
    radii = np.array(radii, ndmin=1, copy=False)

    if seeds.ndim != 2 or seeds.shape[1] != 2:
        raise TypeError('seeds : (nx2)-array')

    if center.ndim != 1 or center.size != 2:
        raise ValueError("center : [x,y]")

    dr2 = (seeds[:, 0] - center[0])**2 + (seeds[:, 1] - center[1])**2

    return seeds[(dr2 + radii**2) <= radius**2, :]
