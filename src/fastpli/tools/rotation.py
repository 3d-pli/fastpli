import numpy as np


def x(phi):
    return np.array(((1, 0, 0), (0, np.cos(phi), -np.sin(phi)),
                     (0, np.sin(phi), np.cos(phi))), np.float32)


def y(phi):
    return np.array(((np.cos(phi), 0, np.sin(phi)), (0, 1, 0),
                     (-np.sin(phi), 0, np.cos(phi))), np.float32)


def z(phi):
    return np.array(((np.cos(phi), -np.sin(phi), 0),
                     (np.sin(phi), np.cos(phi), 0), (0, 0, 1)), np.float32)


def z_2d(phi):
    return z(phi)[:-1, :-1]


def zyz(alpha, beta, gamma):
    return np.dot(z(alpha), np.dot(y(beta), z(gamma)))


def phi(phi):
    return z(phi)


def theta(theta):
    return y(theta)


def theta_phi(theta_, phi_):
    return np.dot(phi(phi_), np.dot(theta(theta_), phi(-phi_)))


def euler(psi, theta_, phi_):
    return np.dot(phi(psi), np.dot(theta(theta_), phi(phi_)))


def a_on_b(a, b):
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
