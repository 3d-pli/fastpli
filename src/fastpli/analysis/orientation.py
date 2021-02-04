# -*- coding: utf-8 -*-
"""
Analyse Methods fiber_bundles orientations
"""

import numpy as np
import numba


@numba.njit(cache=True)
def _remap_direction(phi):
    phi = phi % np.pi
    phi[phi < 0] += np.pi
    return phi


def remap_direction(phi):
    phi = np.array(phi, copy=False)
    shape = phi.shape
    phi.shape = (-1)
    phi = _remap_direction(phi)
    phi.shape = shape
    return phi


@numba.njit(cache=True)
def _remap_orientation(phi, theta):
    phi = phi % (2 * np.pi)
    theta = theta % np.pi

    phi[phi < 0] += 2 * np.pi

    phi[theta < 0] += np.pi
    theta = np.abs(theta)

    mask = theta > 0.5 * np.pi
    phi[mask] += np.pi
    theta[mask] = np.pi - theta[mask]

    phi = phi % (2 * np.pi)

    return phi, theta


def remap_orientation(phi, theta):
    phi = np.array(phi)
    theta = np.array(theta)
    shape = phi.shape

    phi.shape = (-1)
    theta.shape = (-1)

    phi, theta = _remap_orientation(phi, theta)

    phi.shape = shape
    theta.shape = shape

    return phi, theta


@numba.njit(cache=True)
def _remap_spherical(phi, theta):
    phi = phi % (2 * np.pi)
    theta = theta % (2 * np.pi)

    phi[phi < 0] += 2 * np.pi

    phi[theta < 0] += np.pi
    theta = np.abs(theta)

    mask = theta > np.pi
    phi[mask] += np.pi
    theta[mask] = 2 * np.pi - theta[mask]

    phi = phi % (2 * np.pi)

    return phi, theta


def remap_spherical(phi, theta):
    phi = np.array(phi)
    theta = np.array(theta)
    shape = phi.shape

    phi.shape = (-1)
    theta.shape = (-1)

    phi, theta = _remap_spherical(phi, theta)

    phi.shape = shape
    theta.shape = shape

    return phi, theta


def fiber_bundles(fiber_bundles):
    """
    Calculates the orientation of all fiber segments and plots the result.

    Parameters
    ----------
    fiber_bundles : [[(nx4)-array]]
        fiber_bundle object

    Returns
    -------
    res : 1d-array, 1d-array
        arrays of spherical coordinates phi and theta for all fiber segments
    """

    gphi = []
    gtheta = []

    for fb in fiber_bundles:
        for f in fb:
            if f.shape[0] <= 1:
                continue

            df = f[1:, :] - f[0:-1, :]
            phi = np.arctan2(df[:, 1], df[:, 0])
            theta = np.arccos(df[:, 2] / np.linalg.norm(df, axis=1))
            gphi.extend(phi)
            gtheta.extend(theta)

    gphi, gtheta = remap_orientation(gphi, gtheta)

    return gphi, gtheta


def histogram(phi,
              theta,
              ax=None,
              n_phi=100,
              n_theta=50,
              weight_area=False,
              fun=lambda x: x,
              cmap="viridis"):
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
    >>> _, ax = plt.subplots(subplot_kw=dict(projection="polar"))
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
    >>> _, ax = plt.subplots(subplot_kw=dict(projection="polar"))
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

    phi, theta = remap_orientation(phi.ravel(), theta.ravel())

    x = np.linspace(0, 2 * np.pi, n_phi + 1, endpoint=True)
    y = np.linspace(0, np.pi / 2, n_theta + 1, endpoint=True)

    # calculate histogram
    hist, _, _ = np.histogram2d(phi, theta, bins=(x, y))

    if weight_area:
        weights = (np.cos(y[:-1]) - np.cos(y[1:])) * (x[1] - x[0])
        hist /= hist.sum()
        for h in hist:
            h[:] = np.multiply(h[:], 1 / weights)

    if ax:
        X, Y = np.meshgrid(x, np.rad2deg(y))
        pc = ax.pcolormesh(X, Y, fun(hist.T), cmap=cmap)
    else:
        pc = None

    return hist, x, y, pc
