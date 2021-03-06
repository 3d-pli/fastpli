#  -*- coding: utf-8 -*-
"""
Author: Daniel Schmitz INM-1, Forschungszentrum Jülich:
https://doi.org/10.3389/fnana.2018.00075
"""

import itertools
import numpy as np
from scipy import optimize

import numba


@numba.njit(cache=True)
def _invertmatrix(m):
    """
    Calculates the inverse of a 3x3 matrix analytically
    """

    #  calculate determinant
    det = m[0, 0] * m[1, 1] * m[2, 2] + m[1, 0] * m[2, 1] * m[0, 2] + m[
        2, 0] * m[1, 2] - m[0, 0] * m[2, 1] * m[1, 2] - m[2, 0] * m[1, 1] * m[
            0, 2] - m[1, 0] * m[0, 1] * m[2, 2]

    #  initialize inverted matrix
    m_inv = np.empty((3, 3))

    m_inv[0, 0] = m[1, 1] * m[2, 2] - m[1, 2] * m[2, 1]
    m_inv[0, 1] = m[0, 2] * m[2, 1] - m[0, 1] * m[2, 2]
    m_inv[0, 2] = m[0, 1] * m[1, 2] - m[0, 2] * m[1, 1]
    m_inv[1, 0] = m[1, 2] * m[2, 0] - m[1, 0] * m[2, 2]
    m_inv[1, 1] = m[0, 0] * m[2, 2] - m[0, 2] * m[2, 0]
    m_inv[1, 2] = m[0, 2] * m[1, 0] - m[0, 0] * m[1, 2]
    m_inv[2, 0] = m[1, 0] * m[2, 1] - m[1, 1] * m[2, 0]
    m_inv[2, 1] = m[0, 1] * m[2, 0] - m[0, 0] * m[2, 1]
    m_inv[2, 2] = m[0, 0] * m[1, 1] - m[0, 1] * m[1, 0]

    m_inv /= det

    return m_inv


def _calc_centered_vals(stack_2d, gain):

    _, num_rots = stack_2d.shape
    centered_stack = np.empty_like(stack_2d)
    trans = np.mean(stack_2d, 1)
    trans_extended = np.tile(trans, (num_rots, 1)).T
    centered_stack = stack_2d / trans_extended - np.ones_like(stack_2d)
    sigma_stack = np.sqrt(gain * (stack_2d / trans_extended**2 + stack_2d**2 /
                                  (trans_extended**3 * num_rots)))

    return centered_stack, 1 / sigma_stack


@numba.njit(cache=True)
def _symmetrize_angles(direction, inclination):
    """
    Warps angles back into PLI coordinate system in case they are out of bonds
    """
    inclination = (
        ((inclination + np.pi / 2.) % np.pi) -
        np.pi / 2.) * np.sign(-(np.floor(direction / np.pi) % 2) + 0.5)
    direction = direction % np.pi
    return direction, inclination


def _create_rotation_matrix(psi_t, tau):
    """
    Creates a rotation matrix into tilting direction psi_t by tilting angle tau

    Parameters
    ----------
    psi_t : float
    tau : float

    Returns
    -------
    res : (3x3)-array
        rotation matrix
    """

    R = np.array([[
        np.cos(tau) * np.cos(psi_t)**2 + np.sin(psi_t)**2,
        (np.cos(tau) - 1) * np.sin(psi_t) * np.cos(psi_t),
        np.cos(psi_t) * np.sin(tau)
    ],
                  [(np.cos(tau) - 1) * np.sin(psi_t) * np.cos(psi_t),
                   np.cos(tau) * np.sin(psi_t)**2 + np.cos(psi_t)**2,
                   np.sin(psi_t) * np.sin(tau)],
                  [
                      -np.cos(psi_t) * np.sin(tau),
                      -np.sin(psi_t) * np.sin(tau),
                      np.cos(tau)
                  ]])
    return R


@numba.njit(cache=True)
def _create_orientation_vector(phi, alpha):
    """
    Returns orientation vector out of direction angle phi and inclination angle
    alpha

    Parameters
    ----------
    phi : float
    alpha : float

    Returns
    -------
    res : (3x1)-array
        orientation vector
    """

    return np.array([
        np.cos(alpha) * np.cos(phi),
        np.cos(alpha) * np.sin(phi),
        np.sin(alpha)
    ])


def _list_all_tilt_matrices(number_of_tilts, tau):
    """
    Returns an array of rotation matrices

    Parameters
    ----------
    number_of_tilts : int
    tau : float

    Returns
    -------
    res : (number_of_tilts, 3, 3)-array
    """

    tiltangles = np.linspace(0, 2 * np.pi, number_of_tilts + 1)[:-1]
    array_of_matrices = np.zeros((number_of_tilts, 3, 3))
    for tilt in range(number_of_tilts):
        array_of_matrices[tilt, :, :] = _create_rotation_matrix(
            tiltangles[tilt], tau)
    return array_of_matrices


@numba.njit(cache=True)
def _get_direction_inclination(vector):
    """
    Returns direction nand inclination angle from orientation vector

    Parameters
    ----------
    vector : (3, 1)-array

    Returns
    -------
    res : (phi, alpha)
    """

    phi = np.arctan(vector[1] / float(vector[0]))
    alpha = np.arcsin(vector[2])
    return phi, alpha


@numba.njit(cache=True)
def _calculate_rotated_fiber_params(phi, alpha, t_rel, array_matrices, tau):
    """
    Returns rotated fiber parameters

    Parameters
    ----------
    phi : float
    alpha : float
    t_rel : float
    array_matrices : (3, 3)-array
    tau : float

    Returns
    -------
    res : (rotated directions, rotated inclinations, tilted t_rel)
    """

    number_of_tilts = array_matrices.shape[0]
    orientation_vec = _create_orientation_vector(phi, alpha)
    tiltangles = np.linspace(0, 2 * np.pi, number_of_tilts + 1)[:-1]
    tilted_params = np.zeros((number_of_tilts + 1, 2))
    for tilt in range(number_of_tilts):
        tilted_params[tilt + 1, :] = _get_direction_inclination(
            np.dot(array_matrices[tilt], orientation_vec))
    tilted_params[0, :] = [phi, alpha]
    t_rel_array = np.empty(tiltangles.size + 1)
    t_rel_array[1:] = np.full_like(tiltangles, t_rel / np.cos(tau))
    t_rel_array[0] = t_rel
    return tilted_params.T[0, :], tilted_params.T[1, :], t_rel_array


@numba.njit(cache=True)
def _calc_Int_single_fiber_fitted(phi, alpha, t_rel, num_rotations, dir_offset):
    """
    Returns intensity curves in all tilting directions

    Parameters
    ----------
    phi : float
    alpha : float
    t_rel : float
    num_rotations : int

    Returns
    -------
    res : (len(phi)*num_rotations)-array
    """

    phi, alpha = _symmetrize_angles(phi, alpha)
    t_rel = np.abs(t_rel)
    number_tilts = len(phi)
    rotation_angles = np.linspace(0, np.pi, num_rotations + 1)[:-1] + dir_offset

    intensity = np.empty((number_tilts, num_rotations))
    for j in range(0, number_tilts):
        intensity[j, :] = np.sin(np.pi / 2 * t_rel[j] * np.cos(alpha[j])**
                                 2) * np.sin(2 * (rotation_angles - phi[j]))
    return intensity.flatten()


@numba.njit(cache=True)
def _calc_Jacobi(phi, alpha, t_rel, num_rotations, dir_offset):
    """
    Returns Jacobi matrix of the intensity function in all tilting directions

    Parameters
    ----------
    phi : float
    alpha : float
    t_rel : float
    num_rotations : int
    dir_offset : float

    Returns
    -------
    res : (3, len(phi)*num_rotations)-array
    """

    phi, alpha = _symmetrize_angles(phi, alpha)
    t_rel = np.abs(t_rel)
    number_tilts = len(phi)
    rotation_angles = np.linspace(0, np.pi, num_rotations + 1)[:-1] + dir_offset

    Dphi = np.empty((number_tilts, num_rotations))
    Dalpha = np.empty((number_tilts, num_rotations))
    Dtrel = np.empty((number_tilts, num_rotations))

    # calculate derivatives
    for j in range(0, number_tilts):
        Dphi[j, :] = -2 * np.cos(2 * (rotation_angles - phi[j])) * np.sin(
            np.pi / 2 * t_rel[j] * np.cos(alpha[j])**2)
        Dalpha[j, :] = -np.pi * t_rel[j] * np.cos(alpha[j]) * np.sin(
            alpha[j]) * np.cos(np.pi / 2 * t_rel[j] * np.cos(alpha[j])**
                               2) * np.sin(2 * (rotation_angles - phi[j]))
        Dtrel[j, :] = np.pi / 2 * np.cos(alpha[j])**2 * np.sin(
            2 * (rotation_angles - phi[j])) * np.cos(
                np.pi / 2 * t_rel[j] * np.cos(alpha[j])**2)

    # put Jacobian together
    Jacobi = np.empty((3, number_tilts * num_rotations))
    Jacobi[0, :] = Dphi.flatten()
    Jacobi[1, :] = Dalpha.flatten()
    Jacobi[2, :] = Dtrel.flatten()

    return Jacobi


@numba.njit(cache=True)
def _calc_Jacobi_for_opt(phi, alpha, t_rel, num_rotations):
    """
    Returns Jacobi matrix of the intensity function in all tilting directions
    for opt

    Parameters
    ----------
    phi : float
    alpha : float
    t_rel : float
    num_rotations : int

    Returns
    -------
    res : (3, len(phi)*num_rotations)-array
    """

    number_tilts = len(phi)
    rotation_angles = np.linspace(0, np.pi, num_rotations + 1)[:-1]

    Dphi = np.empty((number_tilts, num_rotations))
    Dalpha = np.empty((number_tilts, num_rotations))
    Dtrel = np.empty((number_tilts, num_rotations))

    # calculate derivatives
    for j in range(0, number_tilts):
        Dphi[j, :] = -2 * np.cos(2 * (rotation_angles - phi[j])) * np.sin(
            np.pi / 2 * t_rel[j] * np.cos(alpha[j])**2)
        Dalpha[j, :] = -np.pi * t_rel[j] * np.cos(alpha[j]) * np.sin(
            alpha[j]) * np.cos(np.pi / 2 * t_rel[j] * np.cos(alpha[j])**
                               2) * np.sin(2 * (rotation_angles - phi[j]))
        Dtrel[j, :] = np.pi / 2 * np.cos(alpha[j])**2 * np.sin(
            2 * (rotation_angles - phi[j])) * np.cos(
                np.pi / 2 * t_rel[j] * np.cos(alpha[j])**2)

    # put Jacobian together
    J = np.empty((3,))

    J[0] = np.sum(Dphi)
    J[1] = np.sum(Dalpha)
    J[2] = np.sum(Dtrel)

    return J


def _brute_force_grid(phi, number_inclination_steps, number_t_rel_steps,
                      number_rotations, number_tilt_steps, array_of_matrices,
                      dir_offset, tau, measures, weights):
    """
    Returns the minimumm of the fitted intensities on the brute force grid

    Params:
    float phi - flat direction angle in rad
    int number_inclination_steps - number of inclination steps for gridpoints
    int number_t_reln_steps - number of trel steps for gridpoints
    int number_rotations - number of rotation angles of the measurement
    int number_tilt_steps - number of all measurement poisitions
                            (tilting + 1 flat!)
    float array array array_of_matrices - array of rotation matrices for this
                                          number of tilting positions
    float tau - tilting angle
    float array measures - centered measurement data
    """

    # calculate number of measurements
    number_measurements = number_rotations * number_tilt_steps

    # initialize grid
    incl_steps = np.linspace(-np.pi / 2, np.pi / 2,
                             number_inclination_steps + 2)[1:-1]
    t_rel_steps = np.linspace(0, 0.9, number_t_rel_steps + 2)[1:-1]

    # initialize output array
    number_gridpoints = number_inclination_steps * number_t_rel_steps

    grid_intensities_output = np.empty((number_gridpoints, number_measurements))

    j = 0

    # initialize list of all incl/trel combinations
    gridlist = list(itertools.product(incl_steps, t_rel_steps))

    # loop over all grid points
    for (incl, trel) in itertools.product(incl_steps, t_rel_steps):

        # calculate tilted angles and trel
        phi_array, alpha_array, t_rel_array = _calculate_rotated_fiber_params(
            phi, incl, trel, array_of_matrices, tau)

        # evaluate function and write into output array
        grid_intensities_output[j, :] = _calc_Int_single_fiber_fitted(
            phi_array, alpha_array, t_rel_array, number_rotations, dir_offset)
        j += 1

    # broadcast measurement array number_gridpoints times so it can be
    # substracted from all grid point evaluations
    measures_repeat = np.broadcast_to(measures,
                                      (number_gridpoints, number_measurements))

    # calculate squared deviation
    square_diff = np.square(weights *
                            (grid_intensities_output - measures_repeat))

    # find index of minimal value
    min_arg = np.argmin(np.sum(square_diff, 1))

    return phi, gridlist[min_arg][0], gridlist[min_arg][1]


@numba.njit(cache=True)
def _model(p, array_of_matrices, tau, number_rotations, dir_offset):
    dir_array, alpha_array, t_rel_array = _calculate_rotated_fiber_params(
        p[0], p[1], p[2], array_of_matrices, tau)
    Intensities = _calc_Int_single_fiber_fitted(dir_array, alpha_array,
                                                t_rel_array, number_rotations,
                                                dir_offset)
    return Intensities.flatten()


@numba.njit(cache=True)
def _weighted_jacobi(p, array_of_matrices, tau, number_rotations, dir_offset, _,
                     weights):
    dir_array, alpha_array, t_rel_array = _calculate_rotated_fiber_params(
        p[0], p[1], p[2], array_of_matrices, tau)
    weights = weights.flatten()
    Jacobi = _calc_Jacobi(dir_array, alpha_array, t_rel_array, number_rotations,
                          dir_offset)

    for _ in range(3):
        Jacobi[_, :] *= weights
    return Jacobi


@numba.njit(cache=True)
def _residuum(p, array_of_matrices, tau, number_rotations, dir_offset, measures,
              weights):
    return weights * (_model(p, array_of_matrices, tau, number_rotations,
                             dir_offset) - measures)


@numba.njit(cache=True)
def _chisq(p, array_of_matrices, tau, number_rotations, dir_offset, measures,
           weights):

    residuals = _residuum(p, array_of_matrices, tau, number_rotations,
                          dir_offset, measures, weights)
    return np.sum(residuals * residuals)


def _execute_fit(phi_start, number_inclination_steps, number_t_rel_steps,
                 dir_offset, tau, measures, gain, grad, global_opt_bool):

    # calculate centered intensities and weights
    centered_I, weights = _calc_centered_vals(measures, gain)

    # extract measurement parameters
    number_tilt_steps, number_rotations = centered_I.shape

    # calculate rotation matrices
    rotation_matrices = _list_all_tilt_matrices(number_tilt_steps - 1, tau)

    # flatten measures and weights
    centered_I = centered_I.flatten()
    weights = weights.flatten()
    arguments = (rotation_matrices, tau, number_rotations, dir_offset,
                 centered_I, weights)

    if global_opt_bool is True:  # global

        result = optimize.differential_evolution(_chisq,
                                                 args=arguments,
                                                 bounds=[(0, np.pi),
                                                         (-0.5 * np.pi,
                                                          0.5 * np.pi), (0, 1)],
                                                 popsize=100)

        final_angles = result.x
        fvalue = result.fun
        niter = result.nit

        std_params = np.ones((3,))

        # symmetrize angles back into PLI coordinate space
        direction_out, inclination_out = _symmetrize_angles(
            final_angles[0], final_angles[1])

        # take absolute value of trel
        t_rel = np.abs(final_angles[2])

    else:  # local

        # find start point via brute force minimization
        # maybe speedup: fix startpoint = (phi, 0, 0.2),
        startpoints = _brute_force_grid(phi_start, number_inclination_steps,
                                        number_t_rel_steps, number_rotations,
                                        number_tilt_steps, rotation_matrices,
                                        dir_offset, tau, centered_I, weights)

        # execute least squares fit
        if grad is True:
            final_angles, cov, info, _, _ = optimize.leastsq(
                _residuum,
                startpoints,
                args=arguments,
                full_output=1,
                maxfev=150,
                Dfun=_weighted_jacobi,
                col_deriv=1)
        else:
            final_angles, cov, info, _, _ = optimize.leastsq(_residuum,
                                                             startpoints,
                                                             args=arguments,
                                                             full_output=1,
                                                             maxfev=150)

        # extract value of the cost function
        fvalue = np.sum(np.square(info['fvec']))

        # extract number of iterations of the optimizer
        niter = info['nfev']

        # calculate confidence intervals
        # symmetrize angles back into PLI coordinate space
        direction_out, inclination_out = _symmetrize_angles(
            final_angles[0], final_angles[1])

        # take absolute value of trel
        t_rel = np.abs(final_angles[2])

        #  calculate reduced chi sqaure
        sigma_sq = fvalue / float((len(centered_I) - 3))

        # calculate covariance matrix
        # evaluate jacobian at found minimum
        dir_array, alpha_array, t_rel_array = _calculate_rotated_fiber_params(
            direction_out, inclination_out, t_rel, rotation_matrices, tau)

        J_min = _calc_Jacobi(dir_array, alpha_array, t_rel_array,
                             number_rotations, dir_offset)

        # tile weights for multiplication with jacobi matrix
        weights_repeat = np.tile(weights, (3, 1))
        weighted_jacobi = weights_repeat * J_min
        cov = _invertmatrix(np.dot(weighted_jacobi,
                                   weighted_jacobi.T)) * sigma_sq

        # errors on fit parameters
        std_params = np.diag(np.abs(cov)**0.5)

    return (direction_out, inclination_out, t_rel), std_params, fvalue, niter
