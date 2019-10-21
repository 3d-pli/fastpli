import numpy as np
import warnings
import scipy.ndimage
import multiprocessing as mp
import ctypes


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
