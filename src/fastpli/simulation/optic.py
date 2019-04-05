import numpy as np

from scipy.misc import imresize
from scipy.ndimage.filters import gaussian_filter


def apply(
        image,
        org_intensity,
        org_pixel_size,  # mu meter
        res_pixel_size,  # mu meter
        delta_sigma=0.71,  # only for LAP!
        gain_factor=3,  # only for LAP!
        cropping=0,
        mask=None,
        resize_mode='F'):

    if not mask:
        mask = np.ones_like(image, bool)

    if res_pixel_size % org_pixel_size != 0:
        print('WARNING: OPTIC: res_pixel_size % org_pixel_size:',
              res_pixel_size, org_pixel_size)

    resize = res_pixel_size / org_pixel_size

    res_image = imresize(
        gaussian_filter(image, delta_sigma * resize),
        1 / resize,
        mode=resize_mode)

    # add noise
    if gain_factor > 0:
        zero_mask = res_image == 0  # data is 0 if light is outside of the tissue due to tilting
        res_image[zero_mask] = np.array(
            np.random.negative_binomial(
                res_image[zero_mask] / (gain_factor - 1), 1 / gain_factor),
            np.float32)
        res_image[zero_mask] = 0

    return res_image


def apply_stack(
        image_stack,
        org_intensity,
        org_pixel_size,  # mu meter
        res_pixel_size,  # mu meter
        delta_sigma=0.71,  # only for LAP!
        gain_factor=3,  # only for LAP!
        cropping=0,  # num pixel 
        mask=None,
        resize_mode='F'):

    image_stack = np.swapaxes(image_stack, 0, 2)

    res_image_stack = np.array([
        apply(img, org_intensity, org_pixel_size, res_pixel_size, delta_sigma,
              gain_factor, cropping, mask, resize_mode) for img in image_stack
    ])

    res_image_stack = np.swapaxes(res_image_stack, 0, 2)

    return res_image_stack
