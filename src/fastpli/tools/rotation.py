# -*- coding: utf-8 -*-
"""
Rotation operations
"""

import numpy as np


def x(phi):
    """
    returns rotation array along x-axis

    Parameters
    ----------
    phi: float
        rotation angle

    Returns
    -------
    res : (3,3)-np.ndarray
        rotation array along x-axis
    """
    return np.array(((1, 0, 0), (0, np.cos(phi), -np.sin(phi)),
                     (0, np.sin(phi), np.cos(phi))), float)


def y(phi):
    """
    returns rotation array along y-axis

    Parameters
    ----------
    phi: float
        rotation angle

    Returns
    -------
    res : (3,3)-np.ndarray
        rotation array along y-axis
    """
    return np.array(((np.cos(phi), 0, np.sin(phi)), (0, 1, 0),
                     (-np.sin(phi), 0, np.cos(phi))), float)


def z(phi):
    """
    returns rotation array along z-axis

    Parameters
    ----------
    phi: float
        rotation angle

    Returns
    -------
    res : (3,3)-np.ndarray
        rotation array along z-axis
    """
    return np.array(((np.cos(phi), -np.sin(phi), 0),
                     (np.sin(phi), np.cos(phi), 0), (0, 0, 1)), float)


def z_2d(phi):
    """
    returns rotation array along z-axis for 2d case

    Parameters
    ----------
    phi: float
        rotation angle

    Returns
    -------
    res : (2,2)-np.ndarray
        rotation array along z-axis
    """
    return np.array(((np.cos(phi), -np.sin(phi)), (np.sin(phi), np.cos(phi))),
                    float)


def zyz(alpha, beta, gamma):
    """
    returns rotation array along z-axis y-axis z-axis

    Parameters
    ----------
    alpha: float
        rotation angle
    beta: float
        rotation angle
    gamma: float
        rotation angle

    Returns
    -------
    res : (3,3)-np.ndarray
        rotation array
    """
    return np.dot(z(alpha), np.dot(y(beta), z(gamma)))


def zymz(theta, phi):
    """
    returns rotation array along z-axis y-axis -z-axis

    Parameters
    ----------
    theta: float
        rotation angle
    phi: float
        rotation angle

    Returns
    -------
    res : (3,3)-np.ndarray
        rotation array
    """
    return np.dot(z(phi), np.dot(y(theta), z(-phi)))


def a_on_b(a, b):
    """
    returns rotation array from vector a on vector b

    Parameters
    ----------
    a: np.ndarray
        vector to rotate
    b: np.ndarray
        vector to rotate on

    Returns
    -------
    res : (3,3)-np.ndarray
        rotation array
    """
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    if np.all(a == b):
        return np.identity(3)

    v = np.cross(a, b)
    s = np.linalg.norm(v)
    c = np.dot(a, b)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    R = np.identity(3, float) + vx + np.dot(vx, vx) * (1 - c) / s**2
    return R
