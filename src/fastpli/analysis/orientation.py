# -*- coding: utf-8 -*-
"""
Analyse Methods fiber_bundles orientations
"""

import numpy as np
from numba import njit
import matplotlib.pyplot as plt
import matplotlib.cm as cm


@njit(cache=True)
def _remap_orientation(phi, theta):
    phi = phi % (2 * np.pi)
    theta = theta % (2 * np.pi)

    phi[phi < 0] += 2 * np.pi
    phi[theta < 0] += np.pi
    theta = np.abs(theta)

    mask = theta != theta % np.pi
    phi[mask] += np.pi
    theta = theta % np.pi

    mask = theta > 0.5 * np.pi
    phi[mask] += np.pi
    theta[mask] = np.pi - theta[mask]

    phi = phi % (2 * np.pi)

    # debug
    # if np.any(phi < 0) or np.any(phi >= 2 * np.pi) or np.any(
    #         theta < 0) or np.any(theta > 0.5 * np.pi):
    #     print(phi[phi < 0])
    #     print(phi[phi >= 2 * np.pi])
    #     print(theta[theta < 0])
    #     print(theta[theta > 0.5 * np.pi])
    #     raise ValueError

    return phi, theta


def remap_orientation(phi, theta):
    phi = np.array(phi, copy=False)
    theta = np.array(theta, copy=False)
    shape = phi.shape

    phi.shape = (-1)
    theta.shape = (-1)

    phi, theta = _remap_orientation(phi, theta)

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
              n_angle=100,
              n_radius=50,
              normed=None,
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
    ax : axes object
        for matplotlib
    n_angle : int
        number of angular segments
    n_radius : int
        number of radii segments
    fun : function 
        function apply to histogram height
    cmap : str
        colormap name

    Returns
    -------
    None

    Examples
    --------
    ::

        _, ax = plt.subplots(subplot_kw=dict(projection="polar"))
        phi = np.random.uniform(0,2*np.pi,100)
        theta = np.random.uniform(0,np.pi,100)
        pc = histogram(phi, theta, ax)[-1]
        plt.colorbar(pc, ax=ax)
        ax.set_rmax(90)
        ax.set_rticks(range(0, 90, 10))
        ax.set_rlabel_position(22.5)
        ax.set_yticklabels([])
        ax.grid(True)
        plt.show()
    """

    phi, theta = remap_orientation(phi.ravel(), theta.ravel())

    abins = np.linspace(0, 2 * np.pi, n_angle)
    rbins = np.linspace(0, 90, n_radius)

    #calculate histogram
    hist, _, _ = np.histogram2d(phi,
                                np.rad2deg(theta),
                                bins=(abins, rbins),
                                normed=normed)

    if ax:
        A, R = np.meshgrid(abins, rbins)
        pc = ax.pcolormesh(A, R, fun(hist.T), cmap="viridis")
    else:
        pc = None

    return hist, np.deg2rad(rbins), abins, pc
