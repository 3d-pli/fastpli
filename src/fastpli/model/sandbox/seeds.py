import numpy as np


def triangular_grid(width,
                    height,
                    spacing,
                    center=False,
                    sort=True,
                    endpoint=True):
    """
    Generated 2d triangular grid of seed points inside [0,width]x[0,height].

    Parameters
    ----------
    width, height : float
        length and height of grid [0,width]x[0,height]
    spacing : float
        distance between seed points
    center : bool, optional
        If false, the seed points will be inside [0,width]x[0,height],
        beginning at (0,0).
        If true, the grid will be inside
        [-width/2,width/2]x[-height/2,height/2] with widths eed point at (0,0).
    sort : bool, optional
        If true, the returning seed points are lexsorted along x,y.
    endpoint : bool, optional
        If false, [0,width)x[0,height) or
        [-width/2,width/2)x[-height/2,height/2).
        If true,  [0,width]x[0,height] or
        [-width/2,width/2]x[-height/2,height/2].

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
        width = width / 2
        height = height / 2

    if endpoint:
        if width % spacing == 0 or (width + dx) % spacing == 0:
            width += dx / 2
        if height % spacing == 0 or (height + dy) % spacing == 0:
            height += dy / 2

    grid_0 = np.mgrid[x0:width:spacing, y0:height:2 * dy].reshape(2, -1)
    grid_1 = np.mgrid[x0 + dx:width:spacing,
                      y0 + dy:height:2 * dy].reshape(2, -1)
    grid = np.concatenate((grid_0, grid_1), axis=1)

    if center:
        # mirror x axis
        grid_mirror = grid[:, grid[1, :] != 0]
        grid_mirror[1, :] *= -1
        grid = np.concatenate((grid, grid_mirror), axis=1)
        # mirror y axis
        grid_mirror = grid[:, grid[0, :] != 0]
        grid_mirror[0, :] *= -1
        grid = np.concatenate((grid, grid_mirror), axis=1)

    if sort:
        idx = np.lexsort((grid[0, :], grid[1, :]))
        grid = grid[:, idx]

    return np.ascontiguousarray(grid.T)


def triangular_circle(radius, spacing, center=(0, 0), radii=0):
    """
    Generated 2d triangular grid inside a circle.

    Parameters
    ----------
    radius : float
        radius of circle
    spacing : float
        distance between seed points
    center : (2,)-array_like
        center of circle
    radii : float or (n,)-array_like, optional
        seed points will be iterpreted as cricles with a global or
        individual radii

    Returns
    -------
    res : (nx2)-array
        seed points
    """

    seeds = triangular_grid(2 * (radius + spacing), 2 *
                            (radius + spacing), spacing, True) + center

    return crop_circle(radius, seeds, center, radii)


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
        seed points will be iterpreted as cricles with a global or
        individual radii

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


def crop_circle(radius, seeds, center=(0, 0), radii=0):
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
        seed points will be iterpreted as cricles with a global or
        individual radii

    Returns
    -------
    res : (nx2)-array
        cropped seed points
    """

    seeds = np.array(seeds, ndmin=2, copy=False)
    radii = np.array(radii, ndmin=1, copy=False)

    if seeds.ndim != 2 or seeds.shape[1] != 2:
        raise TypeError('seeds : (nx2)-array')

    return seeds[(np.sum((seeds - center)**2, 1)) <= (radius - radii.T)**2]
