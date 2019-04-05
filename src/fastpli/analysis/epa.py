import numpy as np


def calc(data):

    N = len(data)
    rho = np.deg2rad(np.linspace(0, 180, N, False))

    a0 = 0.0
    a1 = 0.0
    b1 = 0.0

    for i in range(N):
        a0 += data[i]
    a0 = a0 / len(rho)

    for i in range(N):
        a1 += data[i] * np.sin(2 * rho[i])
    a1 = a1 * 2 / N

    for i in range(N):
        b1 += data[i] * np.cos(2 * rho[i])
    b1 *= 2 / N

    t = 2 * a0
    d = 0.5 * np.arctan2(-b1, a1) + np.pi
    r = np.sqrt(a1 * a1 + b1 * b1) / (a0 + 1e-16)

    if d > np.pi:
        d -= np.pi

    return t, d, r


def calc_direction(data):

    N = len(data)
    rho = np.deg2rad(np.linspace(0, 180, N, False))

    a1 = 0.0
    b1 = 0.0

    for i in range(N):
        a1 += data[i] * np.sin(2 * rho[i])
    a1 = a1 * 2 / N

    for i in range(N):
        b1 += data[i] * np.cos(2 * rho[i])
    b1 *= 2 / N

    d = 0.5 * np.arctan2(-b1, a1) + np.pi

    if d > np.pi:
        d -= np.pi

    return d


def calc_direction_map(data):

    N = data.shape[2]
    rho = np.deg2rad(np.linspace(0, 180, N, False))

    sin_2_rho = np.sin(2 * rho)
    cos_2_rho = np.cos(2 * rho)

    data_sin_r_rho = data * sin_2_rho[None, None, :]
    data_cos_r_rho = data * cos_2_rho[None, None, :]

    a1 = np.sum(data_sin_r_rho, 2) * 2 / N
    b1 = np.sum(data_cos_r_rho, 2) * 2 / N

    d = np.array(0.5 * np.arctan2(-b1, a1) + np.pi, data.dtype)

    d[d > np.pi] -= np.pi

    return d


def calc_direction_map(data):

    N = data.shape[2]
    rho = np.deg2rad(np.linspace(0, 180, N, False))

    sin_2_rho = np.sin(2 * rho)
    cos_2_rho = np.cos(2 * rho)

    data_sin_r_rho = data * sin_2_rho[None, None, :]
    data_cos_r_rho = data * cos_2_rho[None, None, :]

    a1 = np.sum(data_sin_r_rho, 2) * 2 / N
    b1 = np.sum(data_cos_r_rho, 2) * 2 / N

    d = np.array(0.5 * np.arctan2(-b1, a1) + np.pi, data.dtype)

    d[d > np.pi] -= np.pi

    return d
