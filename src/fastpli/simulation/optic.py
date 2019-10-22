import numpy as np
import warnings
import scipy.ndimage
from numba import njit


def add_noise(image, gain):
    mask = image > 0
    image[mask] = np.random.negative_binomial(image[mask] / (gain - 1),
                                              1 / gain)
    return image


def filter(image, sigma):
    return scipy.ndimage.filters.gaussian_filter(image, sigma)


def filter2d(image, sigma):
    return np.fft.ifft2(
        scipy.ndimage.fourier_gaussian(np.fft.fft2(image), sigma=sigma)).real


@njit(cache=True)
def resample(image, scale):
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
