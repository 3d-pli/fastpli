# -*- coding: utf-8 -*-
"""
Methods for calculating and applying affine transformation to coordinates and
images.
"""

import numpy as np
import scipy.interpolate

import numba


def _replace_mat_row(B, r, d):
    return np.linalg.det(np.delete(np.vstack([r, B]), (d + 1), axis=0))


@numba.njit(cache=True)
def _nearest_neighbors(image, M):
    """ written for simpli images[x,y,rho]
    """
    image = np.atleast_3d(image)
    image_nn = np.empty_like(image)
    M = np.ascontiguousarray(np.linalg.inv(M))

    x_max, y_max = image.shape[0] - 1, image.shape[1] - 1
    for i in range(image.shape[0]):
        for j in range(image.shape[1]):
            x, y, _ = M @ np.array([i, j, 1.0])
            ii = max(0, min(x_max, int(np.rint(x))))
            jj = max(0, min(y_max, int(np.rint(y))))
            image_nn[i, j, :] = image[ii, jj, :]
    return image_nn


def _interpolate_griddata(image, M, mode):
    """ written for simpli images[x,y,rho]
    """
    image = np.atleast_3d(image)
    image_nn = np.empty_like(image)

    grid_i, grid_j = np.mgrid[0:image.shape[0], 0:image.shape[1]]

    # points -> coordinates in transformed image
    points = np.array(
        [grid_i.flatten(),
         grid_j.flatten(),
         np.ones(grid_j.size)])
    points = (M @ points)[0:2, :]

    for k in range(image.shape[2]):
        image_nn[:, :, k] = scipy.interpolate.griddata(points.T,
                                                       image[:, :, k].flatten(),
                                                       (grid_i, grid_j),
                                                       method=mode)

    return image_nn


def calc_matrix(p_in, p_out):
    """
    Calculate the affine transformation matrix.

    Parameters
    ----------
    p_in, p_out : (3,2)-array_like
        list of 3 x 2d points which will be transformed from p_in to p_out

    Returns
    -------
    res : (3x3)-array
        affine transformation matrix
    """
    p_in = np.array(p_in)
    p_out = np.array(p_out)

    if not np.all(np.equal(np.array(p_in.shape), np.array(p_out.shape))):
        raise TypeError('in and out not the same shape')

    if not np.all(np.equal(np.array(p_in.shape), np.array([3, 2]))):
        print(p_in.shape)
        raise TypeError('shape error: input required [3x2], [3x2]')

    L = p_in.shape[0]
    B = np.vstack([np.transpose(p_in), np.ones(L)])
    D = 1.0 / np.linalg.det(B)
    M = np.array([[(-1)**i * D * _replace_mat_row(B, R, i)
                   for i in range(L)]
                  for R in np.transpose(p_out)])

    return np.vstack([M, [0, 0, 1]])


def exec_matrix(M, x, y):
    """
    Execute the affine transformation.

    Parameters
    ----------
    M : (3,3)-array
        affine transformation matrix
    x, y : float
        2d coordinates to transform

    Returns
    -------
    res : float, float
        transformed coordinates
    """
    x, y, _ = M @ np.array([x, y, 1.0])
    return x, y


def image(image, M, mode='nearest'):
    """
    Execute the affine transformation on simpli images[x,y,rho].

    Parameters
    ----------
    image : 2d-array
        image to transform
    M : float
        affine transformation matrix
    mode : str
        'nearest', 'linear', 'cubic' interpolation mode

    Returns
    -------
    res : 2d-array
        transformed image
    """
    if mode == 'nearest':
        # this is faster then scipy.interpolate.griddata('nearest')
        new_image = _nearest_neighbors(image, M)
    elif mode == 'linear' or mode == 'cubic':
        new_image = _interpolate_griddata(image, M, mode)

    else:
        raise ValueError(f'mode \'{mode}\' does not exist')

    return np.squeeze(new_image)
