import numpy as np
import warnings
import scipy.ndimage
from numba import njit


def add_noise(image, gain, mask=None):
    """
    Adds noise to a simulated ccd image according to the gain value

    Parameters
    ----------
    image : ndarray
        image data. Noise is only added to image values > 0
    gain : float
        gain value for noise level
    mask : ndarray
        noise is only added to masked area

    Returns
    -------
    res : ndarray
    """

    image = np.array(image)

    if not mask:
        mask = image > 0
    else:
        mask = np.logical_and(mask, image > 0)

    image[mask] = np.random.negative_binomial(image[mask] / (gain - 1),
                                              1 / gain)
    return image


def filter(image, sigma):
    """
    Filters image with a 2d gausian kernel

    Parameters
    ----------
    image : ndarray
        image data. Noise is only added to image values > 0
    sigma : float
        variance of gaus kernel

    Returns
    -------
    res : ndarray
    """

    if sigma == 0:
        return image

    return scipy.ndimage.filters.gaussian_filter(image, sigma)


def filter2d(image, sigma):
    """
    Filters image with a 2d gausian kernel

    Parameters
    ----------
    image : ndarray
        image data. Noise is only added to image values > 0
    sigma : float
        variance of gaus kernel

    Returns
    -------
    res : ndarray

    Hints
    -----
    the 2d filter function can be faster then the 1d function.
    """

    if sigma == 0:
        return image

    return np.fft.ifft2(
        scipy.ndimage.fourier_gaussian(np.fft.fft2(image), sigma=sigma)).real


@njit(cache=True)
def resample(image, scale):
    """
    Resamples image acording to scale.
    Resampling is done according to ccd pixels

    Parameters
    ----------
    image : ndarray
        image data. Noise is only added to image values > 0
    scale : float
        scale value for image transformation

    Returns
    -------
    res : ndarray
    """

    sx = int(round(image.shape[0] * scale))
    sy = int(round(image.shape[1] * scale))
    data = np.zeros((sx, sy), dtype=image.dtype)
    divisor = np.zeros((sx, sy), dtype=np.int64)

    for i in range(image.shape[0]):
        ii = int(np.floor((i + 0.5) * scale))
        if ii >= data.shape[0]:
            continue
        for j in range(image.shape[1]):
            jj = int(np.floor((j + 0.5) * scale))
            if jj >= data.shape[1]:
                continue
            divisor[ii, jj] += 1
            data[ii, jj] += image[i, j]
    data = np.divide(data, divisor)

    return data


def resize(image, scale, order=1):
    """
    Resizes image acording to scale

    Parameters
    ----------
    image : ndarray
        image data. Noise is only added to image values > 0
    sigma : float
        variance of gaus kernel

    Returns
    -------
    res : ndarray
    """

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            r"From scipy 0.13.0, the output shape of zoom\(\) is calculated "
            r"with round\(\) instead of int\(\) - for these inputs the size of "
            r"the returned array has changed.")
        return scipy.ndimage.zoom(image, scale, order=order)


def filter_resize(image, delta_sigma, scale, order):
    return resize(filter(image, delta_sigma / scale), scale, order)


def filter_resample(image, delta_sigma, scale):
    return resample(filter(image, delta_sigma / scale), scale)
