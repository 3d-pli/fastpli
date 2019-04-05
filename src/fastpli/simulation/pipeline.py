import numpy as np

from fastpli.analysis import rofl
from fastpli.simulation import optic


def optic_and_rofl(
        image_stack,
        org_intensity,
        org_pixel_size,  # mu meter
        res_pixel_size,  # mu meter
        delta_sigma=0.71,  # only for LAP!
        gain_factor=3,  # only for LAP!
        cropping=0,  # num pixel 
        mask=None,
        resize_mode='F',
        tilt_angle=5.5  # in degree
):

    res_image_stack = optic.apply_stack(
        image_stack, org_intensity, org_pixel_size, res_pixel_size, delta_sigma,
        gain_factor, cropping, mask, resize_mode)

    direction_map = calc_direction_map(res_image_stack)

    rofl_direction, rofl_incl, rofl_t_rel, _, _, _, _, _ = rofl_map(
        res_image_stack, direction_map, tilt_angle, gain_value)

    return rofl_direction, rofl_incl, rofl_t_rel, res_image_stack
