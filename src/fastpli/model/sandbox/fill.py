import numpy as np


def rectangle(a, b, spacing, mode='center'):

    if mode == 'center':
        x0 = -a / 2.0
        y0 = -b / 2.0
    elif mode == 'origin':
        x0 = 0
        y0 = 0
    else:
        ValueError('mode has to be "center" or "origin"')

    dx = np.sin(np.deg2rad(30)) * spacing
    dy = np.cos(np.deg2rad(30)) * spacing

    points_0 = np.mgrid[x0:x0 + a + spacing / 2:spacing, y0:y0 + b +
                        spacing / 2:spacing].reshape(2, -1)
    points_1 = np.mgrid[x0 + dx:x0 + spacing / 2 + a:spacing, y0 + dy:y0 + b +
                        spacing / 2:spacing].reshape(2, -1)

    points = np.concatenate((points_0, points_1), axis=1).T
    points = np.concatenate((points, np.zeros((points.shape[0], 1))), axis=1)

    return points


def circle(radius, spacing):

    dx = np.sin(np.deg2rad(30)) * spacing
    dy = np.cos(np.deg2rad(30)) * spacing

    points_0 = np.mgrid[-radius:radius + spacing / 2:spacing, -radius:radius +
                        spacing / 2:spacing].reshape(2, -1)
    points_1 = np.mgrid[-radius + dx:radius + spacing / 2:spacing, -radius +
                        dy:radius + spacing / 2:spacing].reshape(2, -1)

    points = np.concatenate((points_0, points_1), axis=1)
    dr2 = points[0, :]**2 + points[1, :]**2

    points = points[:, dr2 <= radius**2].T
    points = np.concatenate((points, np.zeros((points.shape[0], 1))), axis=1)

    return points
