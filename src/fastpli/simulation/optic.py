import numpy as np
import warnings
import scipy.ndimage


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


def apply(
        image_stack,
        org_pixel_size,  # mu meter
        res_pixel_size,  # mu meter
        delta_sigma=0.71,  # only for LAP!
        gain=3,  # only for LAP!
        order=1,
        num_threads=2):

    image_stack = np.atleast_3d(np.array(image_stack))
    if image_stack.ndim > 3:
        raise TypeError("image_stack can be 1d, 2d or 3d")

    scale = org_pixel_size / res_pixel_size
    size = np.array(np.round(np.array(image_stack.shape[0:2]) * scale),
                    dtype=int)

    res_image_stack = np.empty((size[0], size[1], image_stack.shape[2]),
                               dtype=image_stack.dtype)

    for i in range(image_stack.shape[2]):
        res_image_stack[:, :, i] = resize(
            filter(image_stack[:, :, i], delta_sigma / scale), scale, order)

    if np.min(res_image_stack.flatten()) < 0:
        raise ValueError("intensity < 0 detected")

    # add noise
    if gain > 0:
        res_image_stack = add_noise(res_image_stack, gain)

    return np.squeeze(res_image_stack)
