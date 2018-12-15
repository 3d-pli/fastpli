import numpy as np


def rot_2d(psi):
    return np.array([[np.cos(psi), -np.sin(psi)],
                     [np.sin(psi), np.cos(psi)]], np.float32)


def rot_phi(psi):
    return np.array([[np.cos(psi), -np.sin(psi), 0],
                     [np.sin(psi), np.cos(psi), 0], [0, 0, 1]], np.float32)


def rot_theta(psi):
    return np.array([[np.cos(psi), 0, np.sin(psi)], [0, 1, 0],
                     [-np.sin(psi), 0, np.cos(psi)]], np.float32)


def rot_theta_phi(theta, phi):
    return np.dot(rot_phi(phi), np.dot(rot_theta(theta), rot_phi(-phi)))


def rot_euler(psi, theta, phi):
    return np.dot(rot_phi(psi), np.dot(rot_theta(theta), rot_phi(phi)))
