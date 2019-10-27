import numpy as np


def triangular_grid(a, b, spacing, center=False, sort=True):
    """
    Generated 2d triangular grid of seed points inside [0,a)x[0,b).

    Parameters
    ----------
    a, b : float
        length and height of grid [0,a)x[0,b)
    spacing : float
        distance between seed points
    center : bool, optional
        If false, the seed points will be inside [0,a)x[0,b), beginning at (0,0)
        If true, the grid will be inside [-a/2,a/2)x[-b/2,b/2) with a seed point at (0,0) 
    sort : bool, optional
        If true, the returning seed points are lexsorted along x,y

    Returns
    -------
    res : (nx2)-array
        seed points
    """

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
    """
    Crops a sequence of 2-dim points inside a rectangle.

    Parameters
    ----------
    a, b : float
        cropping between [0,a]x[0,b]
    a, b : (2,)-array like
        cropping between [a[0],b[0]]x[a[1],b[1]]
    seeds : (n,2)-array like
        to be cropped seed points
    radii : float or (n,)-array like, optional
        seed points will be iterpreted as cricles with a ... or individual radii       
    
    Returns
    -------
    res : (nx2)-array
        cropped seed points
    """

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
    """
    Crops a sequence of 2-dim points inside a circle.

    Parameters
    ----------
    radius : float
        radius of circle area
    seeds : (n,2)-array like
        to be cropped seed points
    center : (2,)-array lile
        center of circle
    radii : float or (n,)-array like, optional
        seed points will be iterpreted as cricles with a ... or individual radii       
    
    Returns
    -------
    res : (nx2)-array
        cropped seed points
    """

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
