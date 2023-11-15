import typing as typ

import matplotlib.pyplot as plt
import numba
import numpy as np

from . import _sh


@numba.njit(cache=True)
def get_num_coeff(bands: int) -> int:
    return ((bands // 2) + 1) * (2 * (bands // 2) + 1)


@numba.njit(cache=True)
def _get_bands_from_coeff(coeff: int) -> int:
    r = 1 / 2 * (np.sqrt(8 * coeff + 1) - 3)
    return int(r)


def _analytic_single_odf(
    cos_theta: np.ndarray, sin_theta: np.ndarray, phi: np.ndarray, bands: int
) -> np.ndarray:
    """original: DOI:10.1109/ISBI.2018.8363804
    This is a helper function for the visualization.
    Currently the visualization does not want to use the _analytic_odf function
    """
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
def _compute_coefficients(
    direction: np.ndarray, inclination: np.ndarray, mask: np.ndarray, bands: int
) -> np.ndarray:
    """Compute odf coefficients from all values in the input arrays

    Args:
        direction (np.ndarray): direction in rad
        inclination (np.ndarray): inclination in rad
        mask (np.ndarray): boolian array of valid entries
        bands (int): _description_

    Returns:
        np.ndarray: _description_
    """

    assert direction.ndim == inclination.ndim == mask.ndim == 1
    assert direction.size == inclination.size == mask.size

    direction = direction.ravel()
    inclination = inclination.ravel()
    mask = mask.ravel()

    direction = direction[mask]
    inclination = inclination[mask]

    if direction.size == 0:
        return np.zeros(get_num_coeff(bands), dtype=np.float32)

    theta = np.pi / 2 - inclination

    return _analytic_odf(np.cos(theta), np.sin(theta), direction, bands)


@numba.njit(cache=True)
def _compute_flatten_array(
    direction: np.ndarray, inclination: np.ndarray, mask: np.ndarray, bands: int
) -> np.array:
    """Compute odf coefficients from all values in the input arrays

    Args:
        direction (np.ndarray): direction in rad
        inclination (np.ndarray): inclination in rad
        mask (np.ndarray): boolian array of valid entries
        bands (int): _description_

    Returns:
        np.ndarray: _description_
    """
    coeff = np.empty((direction.shape[0], get_num_coeff(bands)), dtype=np.float32)
    for i in range(0, coeff.shape[0]):
        coeff[i, :] = _compute_coefficients(
            direction[i, :], inclination[i, :], mask[i, :], bands
        )

    return coeff


def compute(
    direction: np.ndarray,
    inclination: np.ndarray,
    mask: np.ndarray | None = None,
    bands: int = 6,
) -> np.ndarray:
    """calculate odf coefficients
    Odf coefficients are calculated by the analytic solution of the spherical harmonics.
    The input data is interpreted, that the last axis is used for calculating the coefficients.
    The return array is therefore of shape (direction.shape[:-1], get_num_coeff(bands)).

    Args:
        direction (np.ndarray): nd-array of directions in radiant.
        inclination (np.ndarray): nd-array of inclinations in radiant.
        mask (np.ndarray | None, optional): nd-array of bool entries of valid entries. Defaults to None.
        bands (int, optional): number of bands for odf coefficient calculation. Defaults to 6.

    Returns:
        np.ndarray: nd-array of odf coefficients
    """

    ndim = direction.ndim
    direction = np.array(direction, copy=False)
    inclination = np.array(inclination, copy=False)
    assert direction.shape == inclination.shape

    if mask is None:
        mask = np.ones_like(direction, dtype=bool)
    else:
        mask = np.array(mask, copy=False, dtype=bool)
        assert mask.shape == direction.shape

    direction_ = direction.reshape((-1, direction.shape[-1]))
    inclination_ = inclination.reshape((-1, inclination.shape[-1]))
    mask_ = mask.reshape((-1, mask.shape[-1]))

    results = _compute_flatten_array(direction_, inclination_, mask_, bands)

    if ndim == 1:
        results = np.squeeze(results, axis=0)
    else:
        results = results.reshape(direction.shape[:-1] + (get_num_coeff(bands),))

    return results


def _flip_y_coeff(n_elm):
    band = 2 * (1 / 4 * (np.sqrt(8 * n_elm + 1) - 3)) + 1

    if band != int(band):
        raise ValueError("input data misses odf indices")

    idx = 0
    neg = []
    for s in range(0, int(band), 2):
        for si in range(1 + 2 * s):
            if si < ((1 + 2 * s) / 2) - 1:
                neg.append(idx)
            idx += 1

    return neg


def coefficients_to_mrtrix_nii(path: str, coefficients: np.ndarray):
    """save odf coefficients to nii to be able to read in mrtrix"""
    import nibabel as nib

    assert coefficients.ndim == 4  # (x, y, z, coefficients)

    # flip y coefficients
    neg = _flip_y_coeff(coefficients.shape[-1])
    coefficients[:, :, :, neg] *= -1

    nib.save(
        nib.Nifti1Image(coefficients, np.identity(4)),
        path,
    )


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


def visualize_odf(coefficients, n_phi, n_theta, scale=1, fig=None, ax=None):
    phi, theta = np.mgrid[0 : 2 * np.pi : n_phi * 1j, 0 : np.pi : n_theta * 1j]

    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    radius = np.empty_like(phi)
    bands = _get_bands_from_coeff(coefficients.size)

    for i, (p, ct, st) in enumerate(
        zip(phi.ravel(), cos_theta.ravel(), sin_theta.ravel())
    ):
        radius.ravel()[i] = (
            np.sum(
                np.multiply(
                    _analytic_single_odf(ct, st, p, bands),
                    coefficients,
                )
            )
            * scale
        )

    x = radius * np.sin(theta) * np.cos(phi)
    y = radius * np.sin(theta) * np.sin(phi)
    z = radius * np.cos(theta)

    flag = False
    if fig is None:
        fig = plt.figure()
    if ax is None:
        ax = fig.add_subplot(111, projection="3d")
        flag = True
    ax.plot_surface(x, y, z, color="b", alpha=0.5)
    if flag:
        _set_axes_equal(ax)

    return fig, ax
