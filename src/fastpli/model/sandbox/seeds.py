import numpy as np


def triangular_grid(a, b, spacing, center=False, sort=True, endpoint=True):
    """
    Generated 2d triangular grid of seed points inside [0,a]x[0,b].

    Parameters
    ----------
    a, b : float
        length and height of grid [0,a]x[0,b]
    spacing : float
        distance between seed points
    center : bool, optional
        If false, the seed points will be inside [0,a]x[0,b], beginning at (0,0)
        If true, the grid will be inside [-a/2,a/2]x[-b/2,b/2] with a seed point at (0,0) 
    sort : bool, optional
        If true, the returning seed points are lexsorted along x,y
    endpoint : bool, optional
        If false, [0,a)x[0,b) or [-a/2,a/2)x[-b/2,b/2)
        If true,  [0,a]x[0,b] or [-a/2,a/2]x[-b/2,b/2]

    Returns
    -------
    res : (nx2)-array
        seed points
    """

    x0 = 0
    y0 = 0
    dx = spacing / 2  # np.sin(np.deg2rad(30))
    dy = spacing * np.sqrt(3) / 2  # np.cos(np.deg2rad(30))

    if center:
        x0 = -((a / 2) // (2 * dx)) * 2 * dx
        y0 = -((b / 2) // (2 * dy)) * 2 * dy
        a = a / 2
        b = b / 2

    if endpoint:
        if a % spacing == 0 or (a + dx) % spacing == 0:
            a += dx / 2
        if b % spacing == 0 or (b + dy) % spacing == 0:
            b += dy / 2

    grid_0 = np.mgrid[x0:a:spacing, y0:b:2 * dy].reshape(2, -1)
    grid_1 = np.mgrid[x0 + dx:a:spacing, y0 + dy:b:2 * dy].reshape(2, -1)
    grid = np.concatenate((grid_0, grid_1), axis=1)

    if sort:
        idx = np.lexsort((grid[0, :], grid[1, :]))
        grid = grid[:, idx]

    return np.ascontiguousarray(grid.T)


def crop_rectangle(a, b, seeds, radii=0):
    """
    Crops a sequence of 2-dim points inside a rectangle.

    Parameters
    ----------
    a, b : float
        cropping between [0,a]x[0,b]
    a, b : (2,)-array_like
        cropping between [a[0],b[0]]x[a[1],b[1]]
    seeds : (n,2)-array_like
        to be cropped seed points
    radii : float or (n,)-array_like, optional
        seed points will be iterpreted as cricles with a global or individual radii       
    
    Returns
    -------
    res : (nx2)-array
        cropped seed points
    """

    seeds = np.array(seeds, ndmin=2, copy=False)
    radii = np.array(radii, ndmin=2, copy=False)

    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        a, b = [0, 0], [a, b]

    if seeds.ndim != 2 or seeds.shape[1] != 2:
        raise TypeError('seeds : (nx2)-array')

    return seeds[np.logical_and(np.all(seeds - radii.T >= a, 1),
                                np.all(seeds + radii.T <= b, 1))]


def crop_circle(radius, seeds, center=[0, 0], radii=0):
    """
    Crops a sequence of 2-dim points inside a circle.

    Parameters
    ----------
    radius : float
        radius of circle area
    seeds : (n,2)-array_like
        to be cropped seed points
    center : (2,)-array_like
        center of circle
    radii : float or (n,)-array_like, optional
        seed points will be iterpreted as cricles with a global or individual radii       
    
    Returns
    -------
    res : (nx2)-array
        cropped seed points
    """

    seeds = np.array(seeds, ndmin=2, copy=False)
    radii = np.array(radii, ndmin=1, copy=False)

    if seeds.ndim != 2 or seeds.shape[1] != 2:
        raise TypeError('seeds : (nx2)-array')

    return seeds[(np.sum((seeds - center)**2, 1) + radii.T**2) <= radius**2]
