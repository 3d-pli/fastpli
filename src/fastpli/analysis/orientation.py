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

    print(phi.shape)
    print(theta.shape)

    phi = phi % (2 * np.pi)
    theta = theta % np.pi

    phi[phi < 0] += 2 * np.pi

    phi[theta < 0] += np.pi
    theta = np.abs(theta)

    mask = theta > .5 * np.pi
    phi[mask] += np.pi
    theta[mask] = np.pi - theta[mask]

    phi = phi % (2 * np.pi)

    # if np.any(phi < 0) or np.any(phi >= 2 * np.pi) or np.any(
    #         theta < 0) or np.any(theta > 0.5 * np.pi):
    #     raise ValueError

    return phi, theta


def remap_orientation(phi, theta):
    phi = np.array(phi, copy=False)
    theta = np.array(theta, copy=False)
    shape = phi.shape

    phi.shape = (-1)
    theta.shape = (-1)

    _remap_orientation(phi, theta)

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


def histogram(phi, theta, ax, n_angle=100, n_radius=50, cmap="viridis"):
    """
    Plot the Orientation angles in a histogram

    Parameters
    ----------
    phi : array_like
        list of azimuthal angles
    theta : array_like
        list of polar angles
    ax : Axes object

    Returns
    -------
    None

    Examples
    --------
    ::

        _, ax = plt.subplots(subplot_kw=dict(projection="polar"))
        phi = np.random(0,2*np.pi,100)
        theta = np.random(0,np.pi,100)
        pc = histogram(phi, theta, ax)
        plt.colorbar(pc, ax=ax)
        plt.show()
    """

    phi, theta = remap_orientation(phi.ravel(), theta.ravel())

    abins = np.linspace(0, 2 * np.pi, n_angle)
    rbins = np.linspace(0, 90, n_radius)

    #calculate histogram
    hist, _, _ = np.histogram2d(phi, np.rad2deg(theta), bins=(abins, rbins))
    A, R = np.meshgrid(abins, rbins)

    pc = ax.pcolormesh(A, R, hist.T, cmap="viridis")

    return pc

    # ax.set_rmax(90)
    # ax.set_rticks(range(0, 90, 10))
    # ax.set_rlabel_position(22.5)
    # ax.set_yticklabels([])
    # ax.grid(True)
