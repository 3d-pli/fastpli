import typing as typ

import matplotlib.pyplot as plt
import numba
import numpy as np

from . import _sh


@numba.njit(cache=True)
def get_num_coeff(bands: int) -> int:
    return ((bands // 2) + 1) * (2 * (bands // 2) + 1)


def _get_bands_from_coeff(coeff: int) -> int:
    r = 1 / 2 * (np.sqrt(8 * coeff + 1) - 3)
    if r % 1 != 0:
        raise ValueError(f"num coeff={coeff} is invalid")
    return int(r)


def _analytic_single_odf(
    cos_theta: np.ndarray, sin_theta: np.ndarray, phi: np.ndarray, bands: int
) -> np.ndarray:
    """original: DOI:10.1109/ISBI.2018.8363804"""
    n_coeff = get_num_coeff(bands)
    real_sph_harm = np.empty(n_coeff, np.float32)

    i = 0
    for b in np.arange(0, bands + 1, 2):
        for o in np.arange(-b, b + 1, 1):
            real_sph_harm[i] = (
                np.sum(_sh.spherical_harmonics(b, o, cos_theta, sin_theta, phi))
                / phi.size
            )
            i += 1

    return real_sph_harm


@numba.njit(cache=True)
def _analytic_odf(
    cos_theta: np.ndarray, sin_theta: np.ndarray, phi: np.ndarray, bands: int
) -> np.ndarray:
    """original: DOI:10.1109/ISBI.2018.8363804"""
    n_coeff = get_num_coeff(bands)
    real_sph_harm = np.empty(n_coeff, np.float32)

    i = 0
    for b in np.arange(0, bands + 1, 2):
        for o in np.arange(-b, b + 1, 1):
            real_sph_harm[i] = (
                np.sum(_sh.spherical_harmonics(b, o, cos_theta, sin_theta, phi))
                / phi.size
            )
            i += 1

    return real_sph_harm


@numba.njit(cache=True)
def calc_real_sh(
    dir_data: np.ndarray, incl_data: np.ndarray, mask_data: np.ndarray, bands: int
) -> np.ndarray:
    assert dir_data.size == incl_data.size == mask_data.size

    dir_data = dir_data.ravel()
    incl_data = incl_data.ravel()
    mask_data = mask_data.ravel()

    dir_data = dir_data[mask_data]
    incl_data = incl_data[mask_data]

    if dir_data.size == 0:
        return np.zeros(get_num_coeff(bands), dtype=np.float32)

    phi_data = np.radians(dir_data)
    theta_data = np.radians(90 - incl_data)

    return _analytic_odf(np.cos(theta_data), np.sin(theta_data), phi_data, bands)


@numba.njit(cache=True)
def compute(
    direction: np.ndarray,
    inclination: np.ndarray,
    mask: typ.Optional[np.ndarray],
    bands: int,
) -> np.ndarray:
    if mask is None:
        mask = np.ones_like(direction, dtype=np.uint8)

    coeff = np.empty((direction.shape[0], get_num_coeff(bands)), dtype=np.float32)

    for x in range(0, coeff.shape[0]):
        coeff[x, :] = calc_real_sh(
            direction[x, :], inclination[x, :], mask[x, :], bands
        )

    return coeff


def _set_axes_equal(ax):
    """
    Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    """

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


def visualize_odf(coefficients, n_phi, n_theta, fig=None, ax=None, show=True):
    phi, theta = np.mgrid[0 : 2 * np.pi : n_phi * 1j, 0 : np.pi : n_theta * 1j]

    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    radius = np.empty_like(phi)
    bands = _get_bands_from_coeff(coefficients.size)

    for i, (p, ct, st) in enumerate(
        zip(phi.ravel(), cos_theta.ravel(), sin_theta.ravel())
    ):
        radius.ravel()[i] = np.sum(
            np.multiply(
                _analytic_single_odf(ct, st, p, bands),
                coefficients,
            )
        )

    x = radius * np.sin(theta) * np.cos(phi)
    y = radius * np.sin(theta) * np.sin(phi)
    z = radius * np.cos(theta)

    if fig is None:
        fig = plt.figure()
    if ax is None:
        ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(x, y, z, color="b", alpha=0.5)
    if show:
        plt.plot([0, 1], [0, 0], [0, 0], color="r")
        plt.plot([0, 0], [0, 1], [0, 0], color="g")
        plt.plot([0, 0], [0, 0], [0, 1], color="b")
        _set_axes_equal(ax)
        plt.show()

    return fig, ax
