# -*- coding: utf-8 -*-
"""
Analyse Methods fiber_bundles orientations
"""
import numpy as np
import numba

# pylint: disable = unpacking-non-sequence
# pylint: disable = no-value-for-parameter


@numba.vectorize([numba.float32(numba.float32),
                  numba.float64(numba.float64)],
                 nopython=True,
                 cache=True)
def remap_direction(phi):
    """
    Return direction in range of [0, np.pi)

    Returns
    -------
    res : np.ndarray
        direction in [0, np.pi)
    """

    phi_ = phi % np.pi
    if phi_ < 0:
        phi_ += np.pi

    return phi_


@numba.guvectorize(
    [(numba.float32[:], numba.float32[:], numba.float32[:], numba.float32[:]),
     (numba.float64[:], numba.float64[:], numba.float64[:], numba.float64[:])],
    '(),()->(),()',
    nopython=True,
    cache=True)
def remap_sphere(phi, theta, phi_, theta_):
    """
    Return the azimuthal angle in range of [0, 2*np.pi) and the polar angle theta in [0, np.pi]

    Returns
    -------
    res : (np.ndarray, np.ndarray)
        phi in [0, 2*np.pi), theta in [0, np.pi]
    """
    phi_[:] = phi % (2 * np.pi)
    theta_[:] = theta % (2 * np.pi)

    phi_[phi_ < 0] += 2 * np.pi

    phi_[theta_ < 0] += np.pi
    theta_[:] = np.abs(theta_)

    mask = theta_ > np.pi
    phi_[mask] += np.pi
    theta_[mask] = 2 * np.pi - theta_[mask]

    phi_[:] = phi_ % (2 * np.pi)
    phi_[theta_ == 0] = 0
    phi_[theta_ == np.pi] = 0


@numba.guvectorize(
    [(numba.float32[:], numba.float32[:], numba.float32[:], numba.float32[:]),
     (numba.float64[:], numba.float64[:], numba.float64[:], numba.float64[:])],
    '(),()->(),()',
    nopython=True,
    cache=True)
def remap_half_sphere_z(phi, theta, phi_, theta_):
    """
    Return the azimuthal angle in range of [0, 2*np.pi) and the polar angle theta in [0, 0.5*np.pi]

    Returns
    -------
    res : (np.ndarray, np.ndarray)
        aximuth in [0, 2*np.pi), theta in [0, 0.5*np.pi]
    """

    # remap to sphere
    phi_[:] = phi % (2 * np.pi)
    theta_[:] = theta % (2 * np.pi)

    phi_[phi_ < 0] += 2 * np.pi

    phi_[theta_ < 0] += np.pi
    theta_[:] = np.abs(theta_)

    mask = theta_ > np.pi
    phi_[mask] += np.pi
    theta_[mask] = 2 * np.pi - theta_[mask]

    # remap to half sphere -> z
    mask = theta_ > np.pi / 2
    phi_[mask] += np.pi
    theta_[mask] -= np.pi / 2

    phi_[:] = phi_ % (2 * np.pi)
    phi_[theta_ == 0] = 0


@numba.guvectorize(
    [(numba.float32[:], numba.float32[:], numba.float32[:], numba.float32[:]),
     (numba.float64[:], numba.float64[:], numba.float64[:], numba.float64[:])],
    '(),()->(),()',
    nopython=True,
    cache=True)
def remap_half_sphere_x(phi, theta, phi_, theta_):
    """
    Return the azimuthal angle in range of [-np.pi, np.pi) and the
    polar angle in [0, np.pi]

    Returns
    -------
    res : (np.ndarray, np.ndarray)
        aximuth in [-np.pi/2, np.pi/2), theta in [0, np.pi]
    """

    # remap to sphere
    phi_[:] = phi % (2 * np.pi)
    theta_[:] = theta % (2 * np.pi)

    phi_[phi_ < 0] += 2 * np.pi

    phi_[theta_ < 0] += np.pi
    theta_[:] = np.abs(theta_)

    mask = theta_ > np.pi
    phi_[mask] += np.pi
    theta_[mask] = 2 * np.pi - theta_[mask]

    phi_[:] = phi_ % (2 * np.pi)
    phi_[theta_ == 0] = 0
    phi_[theta_ == np.pi] = 0

    # transform to half sphere -> x
    mask = np.logical_and(phi_ >= np.pi / 2, phi_ < 3 * np.pi / 2)
    phi_[mask] -= np.pi
    theta_[mask] = np.pi - theta_[mask]

    mask = phi_ >= 3 * np.pi / 2
    phi_[mask] -= 2 * np.pi

    phi_[theta_ == 0] = 0
    phi_[theta_ == np.pi] = 0
    theta_[theta_ == np.pi] = 0


def fiber(fiber):
    """
    Calculates the orientation of all fiber segments.

    Parameters
    ----------
    fiber : (nx4)-array
        fiber object

    Returns
    -------
    res : 1d-array, 1d-array
        arrays of spherical coordinates phi and theta for all fiber segments
    """

    gphi = []
    gtheta = []

    if fiber.shape[0] <= 1:
        return None, None

    df = fiber[1:, :] - fiber[0:-1, :]
    phi = np.arctan2(df[:, 1], df[:, 0])
    theta = np.arccos(df[:, 2] / np.linalg.norm(df, axis=1))
    gphi.extend(phi)
    gtheta.extend(theta)

    gphi, gtheta = remap_half_sphere_z(gphi, gtheta)

    return gphi, gtheta


def fiber_bundle(fb):
    """
    Calculates the orientation of all fiber segments.

    Parameters
    ----------
    fb : [(nx4)-array]
        fiber_bundle object

    Returns
    -------
    res : 1d-array, 1d-array
        arrays of spherical coordinates phi and theta for all fiber segments
    """

    gphi = []
    gtheta = []

    for f in fb:
        if f.shape[0] <= 1:
            continue

        df = f[1:, :] - f[0:-1, :]
        phi = np.arctan2(df[:, 1], df[:, 0])
        theta = np.arccos(df[:, 2] / np.linalg.norm(df, axis=1))
        gphi.extend(phi)
        gtheta.extend(theta)

    gphi, gtheta = remap_half_sphere_z(gphi, gtheta)

    return gphi, gtheta


def fiber_bundles(fbs):
    """
    Calculates the orientation of all fiber segments.

    Parameters
    ----------
    fbs : [[(nx4)-array]]
        fiber_bundle object

    Returns
    -------
    res : 1d-array, 1d-array
        arrays of spherical coordinates phi and theta for all fiber segments
    """

    gphi = []
    gtheta = []

    for fb in fbs:
        for f in fb:
            if f.shape[0] <= 1:
                continue

            df = f[1:, :] - f[0:-1, :]
            phi = np.arctan2(df[:, 1], df[:, 0])
            theta = np.arccos(df[:, 2] / np.linalg.norm(df, axis=1))
            gphi.extend(phi)
            gtheta.extend(theta)

    gphi, gtheta = remap_half_sphere_z(gphi, gtheta)

    return gphi, gtheta


def histogram(phi,
              theta,
              n_phi=100,
              n_theta=50,
              weight_area=False,
              fun=lambda x: x,
              ax=None,
              cmap='viridis'):
    """
    Plot the Orientation angles in a histogram

    Parameters
    ----------
    phi : array_like
        list of azimuthal angles
    theta : array_like
        list of polar angles
    ax : axes object, optional
        for matplotlib
    n_phi : int, optional
        number of angular segments
    n_theta : int, optional
        number of radii segments
    density : bool, optional
        density argument for numpy histogram
    weight_area : bool, optional
        weighting the density by the histogram bin on a sphere
    fun : function , optional
        function apply to histogram height
    cmap : str, optional
        colormap name

    Returns
    -------
    None

    Examples
    --------
    >>> # counts
    >>> _, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    >>> _, _, _, pc = histogram(phi,
                                theta,
                                ax=ax,
                                n_phi=60,
                                n_theta=30,
                                weight_area=False)
    >>> cbar = plt.colorbar(pc, ax=ax)
    >>> cbar.ax.set_title('#')
    >>> ax.set_rmax(90)
    >>> ax.set_rticks(range(0, 90, 10))
    >>> ax.set_rlabel_position(22.5)
    >>> ax.set_yticklabels([])
    >>> ax.grid(True)

    >>> # density
    >>> _, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    >>> phi = np.random.normal(np.pi / 3, 0.5, 1000)
    >>> theta = np.random.normal(np.deg2rad(45), 0.5, 1000)

    >>> _, _, _, pc = histogram(phi,
                                theta,
                                ax=ax,
                                n_phi=60,
                                n_theta=30,
                                weight_area=True)
    >>> cbar = plt.colorbar(pc, ax=ax)
    >>> cbar.ax.set_title('$P(\\vartheta, \\varphi)$')

    >>> ax.set_rmax(90)
    >>> ax.set_rticks(range(0, 90, 10))
    >>> ax.set_rlabel_position(22.5)
    >>> ax.set_yticklabels([])
    >>> ax.set_yticklabels([])
    >>> ax.grid(True)

    >>> plt.show()
    """

    phi, theta = remap_half_sphere_z(phi.ravel(), theta.ravel())

    x = np.linspace(0, 2 * np.pi, n_phi + 1, endpoint=True)
    y = np.linspace(0, np.pi / 2, n_theta + 1, endpoint=True)

    # calculate histogram
    hist, _, _ = np.histogram2d(phi, theta, bins=(x, y))

    if weight_area:
        weights = (np.cos(y[:-1]) - np.cos(y[1:])) * (x[1] - x[0])
        hist /= hist.sum()
        for h in hist:
            h[:] = np.multiply(h[:], 1 / weights)

    hist = fun(hist)

    if ax:
        X, Y = np.meshgrid(x, np.rad2deg(y))
        pc = ax.pcolormesh(X, Y, hist.T, cmap=cmap)
    else:
        pc = None

    return hist, x, y, pc
