import numpy as np


def _vec_to_rgb(x, y, z):
    l = np.sqrt(x**2 + y**2 + z**2)
    x = int(x / l * 255)
    y = int(y / l * 255)
    z = int(z / l * 255)
    return np.array([np.abs(x), np.abs(y), np.abs(z)])


def _hsv_black_to_rgb_space(h, s, v):
    # images have to be saved in rgb space

    h = (h + 360) % 360

    hi = np.floor(h / 60)
    f = h / 60.0 - hi

    p = v * (1 - s)
    q = v * (1 - s * f)
    t = v * (1 - s * (1 - f))

    if hi == 1:
        rgb = np.array((q, v, p))
    elif hi == 2:
        rgb = np.array((p, v, t))
    elif hi == 3:
        rgb = np.array((p, q, v))
    elif hi == 4:
        rgb = np.array((t, p, v))
    elif hi == 5:
        rgb = np.array((v, p, q))
    else:
        rgb = np.array((v, t, p))

    return np.array(rgb * 255, int)


def _orientation_to_hsv(directionValue, inclinationValue):
    h = 360.0 * np.abs(directionValue) / np.pi
    s = 1.0
    v = 1.0 - (2 * np.abs(inclinationValue) / np.pi)

    return _hsv_black_to_rgb_space(h, s, v)


def hsv_black_sphere(n=128):
    sphere = np.zeros((n, n, 3), dtype=np.uint8)

    for x in range(n):
        xx = (-n // 2 + x) / (n // 2)
        for y in range(n):
            yy = (-n // 2 + y) / (n // 2)
            if xx**2 + yy**2 > 1:
                continue
            zz = np.sqrt(1 - xx**2 - yy**2)

            direction = np.arctan2(yy, xx)
            inclination = np.arcsin(zz)

            if direction < 0:
                direction += np.pi

            sphere[x, y, :] = _orientation_to_hsv(direction, inclination)

    return sphere


def rgb_sphere(n=128):
    sphere = np.zeros((n, n, 3), dtype=np.uint8)

    for x in range(n):
        xx = (n // 2 - x) / (n // 2)
        for y in range(n):
            yy = (n // 2 - y) / (n // 2)
            if xx**2 + yy**2 > 1:
                continue
            zz = np.sqrt(1 - xx**2 - yy**2)

            sphere[x, y, :] = _vec_to_rgb(xx, yy, zz)

    return sphere


def unit_vectors(direction, inclination, mask=None):
    UnitX = np.sin(0.5 * np.pi - inclination) * np.cos(direction)
    UnitY = np.sin(0.5 * np.pi - inclination) * np.sin(direction)
    UnitZ = np.cos(0.5 * np.pi - inclination)

    if mask is not None:
        UnitX[~mask] = 0
        UnitY[~mask] = 0
        UnitZ[~mask] = 0

    return UnitX, UnitY, UnitZ


def fom_hsv_black(direction, inclination, mask=None):
    if mask is None:
        mask = np.ones_like(direction, dtype=np.bool)

    hsv = np.zeros((mask.shape[0], mask.shape[1], 3), np.uint8)
    for x in range(mask.shape[0]):
        for y in range(mask.shape[1]):
            if not mask[x, y]:
                continue

            hsv[x, y, :] = _orientation_to_hsv(direction[x, y],
                                               inclination[x, y])
    return hsv


def fom_rgb(direction, inclination, mask=None):
    if mask is None:
        mask = np.ones_like(direction, dtype=np.bool)

    rgb = np.zeros((mask.shape[0], mask.shape[1], 3), np.uint8)
    for x in range(mask.shape[0]):
        for y in range(mask.shape[1]):
            if not mask[x, y]:
                continue

            rgb[x, y, :] = _vec_to_rgb(
                np.sin(0.5 * np.pi - inclination[x, y]) *
                np.cos(direction[x, y]),
                np.sin(0.5 * np.pi - inclination[x, y]) *
                np.sin(direction[x, y]),
                np.cos(0.5 * np.pi - inclination[x, y]))
    return rgb
