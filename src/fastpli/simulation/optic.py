import numpy as np

from PIL import Image
from scipy.ndimage.filters import gaussian_filter


def apply(
        image,
        org_pixel_size,  # mu meter
        res_pixel_size,  # mu meter
        delta_sigma=0.71,  # only for LAP!
        gain=3,  # only for LAP!
        cropping=0,
        resample_mode=Image.BILINEAR):

    if (res_pixel_size / org_pixel_size) % 1.0 != 0:
        print('WARNING: OPTIC: res_pixel_size % org_pixel_size:',
              res_pixel_size, org_pixel_size)  #

    resize = res_pixel_size / org_pixel_size
    new_size = np.array(np.array(image.shape) // resize, dtype=int)

    res_image = np.array(
        Image.fromarray(gaussian_filter(image, delta_sigma * resize)).resize(
            new_size, resample_mode))

    # add noise
    if gain > 0:
        mask = res_image != 0  # data is 0 if light is outside of the tissue due to tilting

        res_image[mask] = np.array(
            np.random.negative_binomial(res_image[mask] / (gain - 1), 1 / gain),
            np.float32)

    return res_image


def resize_img(image,
               org_pixel_size,
               res_pixel_size,
               resample_mode=Image.BILINEAR):
    if res_pixel_size % org_pixel_size != 0:
        print('WARNING: OPTIC: res_pixel_size % org_pixel_size:',
              res_pixel_size, org_pixel_size)

    resize = res_pixel_size / org_pixel_size
    new_size = np.array(np.array(image.shape) // resize, dtype=int)

    res_image = np.array(Image.fromarray(image).resize(new_size, resample_mode))

    return res_image


def apply_stack(
        image_stack,
        org_pixel_size,  # mu meter
        res_pixel_size,  # mu meter
        delta_sigma=0.71,  # only for LAP!
        gain=3,  # only for LAP!
        cropping=0,  # num pixel
        resize_mode='F'):

    image_stack = np.swapaxes(image_stack, 0, 2)

    res_image_stack = np.array([
        apply(img, org_pixel_size, res_pixel_size, delta_sigma, gain, cropping,
              resize_mode) for img in image_stack
    ])

    res_image_stack = np.swapaxes(res_image_stack, 0, 2)

    return res_image_stack
