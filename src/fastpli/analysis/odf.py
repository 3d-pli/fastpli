import numba
import numpy as np
import typing as typ

from . import _sh


@numba.njit(cache=True)
def get_num_coeff(bands: int) -> int:
    return ((bands // 2) + 1) * (2 * (bands // 2) + 1)


@numba.njit(cache=True)
def _analytic_odf(cos_theta: np.ndarray, sin_theta: np.ndarray, phi: np.ndarray,
                  bands: int) -> np.ndarray:
    """original: DOI:10.1109/ISBI.2018.8363804
    """
    n_coeff = get_num_coeff(bands)
    real_sph_harm = np.empty(n_coeff, np.float32)

    i = 0
    for b in np.arange(0, bands + 1, 2):
        for o in np.arange(-b, b + 1, 1):
            real_sph_harm[i] = np.sum(
                _sh.sph_harm(b, o, cos_theta, sin_theta, phi)) / phi.size
            i += 1

    return real_sph_harm


@numba.njit(cache=True)
def calc_real_sh(dir_data: np.ndarray, incl_data: np.ndarray,
                 mask_data: np.ndarray, bands: int) -> np.ndarray:
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

    return _analytic_odf(np.cos(theta_data), np.sin(theta_data), phi_data,
                         bands)


@numba.njit(cache=True)
def compute(direction: np.ndarray, inclination: np.ndarray,
            mask: typ.Optional[np.ndarray], bands: int) -> np.ndarray:

    if mask is None:
        mask = np.ones_like(direction, dtype=np.uint8)

    coeff = np.empty((direction.shape[0], get_num_coeff(bands)),
                     dtype=np.float32)

    for x in range(0, coeff.shape[0]):
        coeff[x, :] = calc_real_sh(direction[x, :], inclination[x, :],
                                   mask[x, :], bands)

    return coeff
