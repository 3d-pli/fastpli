# -*- coding: utf-8 -*-
"""
Rotation operations
"""

import numpy as np


def x(phi):
    """ 3d rotation around x-axis: float -> (3,3)-array """
    return np.array(((1, 0, 0), (0, np.cos(phi), -np.sin(phi)),
                     (0, np.sin(phi), np.cos(phi))), float)


def y(phi):
    """ 3d rotation around y-axis: float -> (3,3)-array """
    return np.array(((np.cos(phi), 0, np.sin(phi)), (0, 1, 0),
                     (-np.sin(phi), 0, np.cos(phi))), float)


def z(phi):
    """ 3d rotation around z-axis: float -> (3,3)-array """
    return np.array(((np.cos(phi), -np.sin(phi), 0),
                     (np.sin(phi), np.cos(phi), 0), (0, 0, 1)), float)


def z_2d(phi):
    """ 2d rotation around z-axis: float -> (2,2)-array """
    return np.array(((np.cos(phi), -np.sin(phi)), (np.sin(phi), np.cos(phi))),
                    float)


def zyz(alpha, beta, gamma):
    """ 3d rotation around z-, y-, z-axis:
    float, float, float -> (3,3)-array
    """
    return np.dot(z(alpha), np.dot(y(beta), z(gamma)))


def zymz(theta, phi):
    """ 3d rotation around (theta,phi)-axis: float, float -> (3,3)-array """
    return np.dot(z(phi), np.dot(y(theta), z(-phi)))


def a_on_b(a, b):
    """
    return rotation matrix when rotating a on b:
    (3)-array, (3)-array -> (3,3)-array
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
