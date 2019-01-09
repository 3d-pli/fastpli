import numpy as np

# TODO: n*pi cases


def x(phi):
    return np.array(((1, 0, 0), (0, np.cos(phi), -np.sin(phi)), (0, np.sin(phi), np.cos(phi))), np.float32)


def y(phi):
    return np.array(((np.cos(phi), 0, np.sin(phi)), (0, 1, 0), (-np.sin(phi), 0, np.cos(phi))), np.float32)


def z(phi):
    return np.array(((np.cos(phi), -np.sin(phi), 0), (np.sin(phi), np.cos(phi), 0), (0, 0, 1)), np.float32)


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
