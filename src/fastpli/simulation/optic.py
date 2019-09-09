import numpy as np
import warnings

from PIL import Image
from scipy.ndimage.filters import gaussian_filter


def add_noise(image, gain, mask=None):
    image = np.copy(image)

    if mask is None:
        mask = image > 0

    image[mask] = np.array(
        np.random.negative_binomial(image[mask] / (gain - 1), 1 / gain), float)

    return image


def filter(image, sigma):
    return gaussian_filter(image, sigma)


def resize(image, size, resample_mode=Image.BILINEAR):
    return np.array(Image.fromarray(image).resize(size, resample_mode))


def apply(
        image_stack,
        org_pixel_size,  # mu meter
        res_pixel_size,  # mu meter
        delta_sigma=0.71,  # only for LAP!
        gain=3,  # only for LAP!
        resample_mode=Image.BILINEAR):

    image_stack = np.atleast_3d(np.array(image_stack))
    if image_stack.ndim > 3:
        raise TypeError("image_stack can be 1d, 2d or 3d")

    if (res_pixel_size / org_pixel_size) % 1.0 != 0:
        warnings.warn('OPTIC: res_pixel_size % org_pixel_size: ' +
                      str(res_pixel_size) + '%' + str(org_pixel_size))

    scale = res_pixel_size / org_pixel_size
    size = np.array(np.array(image_stack.shape[0:2]) // scale, dtype=int)

    res_image_stack = np.empty((size[0], size[1], image_stack.shape[2]),
                               dtype=image_stack.dtype)

    for i in range(image_stack.shape[2]):
        res_image_stack[:, :, i] = resize(
            filter(image_stack[:, :, i], delta_sigma * scale), size,
            resample_mode)

    # add noise
    if gain > 0:
        for i in range(image_stack.shape[2]):
            res_image_stack[:, :, i] = add_noise(res_image_stack[:, :, i], gain)

    return np.squeeze(res_image_stack)
