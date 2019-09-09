import numpy as np
import warnings

from PIL import Image
from scipy.ndimage.filters import gaussian_filter


def add_noise(image, gain, mask=None):
    image = np.copy(image)

    if mask is None:
        mask = image > 0

    image[mask] = np.array(
        np.random.negative_binomial(res_image[mask] / (gain - 1), 1 / gain),
        float)

    return image


def filter(image, sigma):
    return gaussian_filter(image, sigma)


def resize(image, size, resample_mode=Image.BILINEAR):
    return np.array(Image.fromarray(image).resize(size, resample_mode))


def apply(
        image,
        org_pixel_size,  # mu meter
        res_pixel_size,  # mu meter
        delta_sigma=0.71,  # only for LAP!
        gain=3,  # only for LAP!
        resample_mode=Image.BILINEAR):

    image = np.copy(image)

    if (res_pixel_size / org_pixel_size) % 1.0 != 0:
        warnings.warn('OPTIC: res_pixel_size % org_pixel_size: ' +
                      str(res_pixel_size) + '%' + str(org_pixel_size))

    scale = res_pixel_size / org_pixel_size
    size = np.array(np.array(image.shape) // scale, dtype=int)
    image = resize(filter(image, delta_sigma * scale), size, resample_mode)

    # add noise
    if gain > 0:
        image = add_noise(image, gain)

    return image
