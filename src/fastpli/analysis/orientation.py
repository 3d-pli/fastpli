import numpy as np
from numba import njit


@njit(cache=True)
def _remap_orientation(phi, theta):
    phi = phi % (2 * np.pi)
    theta = theta % np.pi

    phi[phi < 0] += 2 * np.pi

    phi[theta < 0] += np.pi
    theta = np.abs(theta)

    phi[theta > .5 * np.pi] += np.pi
    theta[theta > .5 * np.pi] = np.pi - theta[theta >= .5 * np.pi]

    phi = phi % (2 * np.pi)

    if np.any(phi < 0) or np.any(phi >= 2 * np.pi) or np.any(
            theta < 0) or np.any(theta > 0.5 * np.pi):
        raise ValueError

    return phi, theta


def fiber_statistic(fiber_bundles,
                    plot='scatter',
                    rticks=[10, 20, 30, 40, 50, 60, 70, 80, 90]):

    if plot:
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm

    gphi = []
    gtheta = []

    fig, ax = plt.subplots(subplot_kw=dict(projection="polar"))
    colors = cm.viridis(np.linspace(0, 1, len(fiber_bundles)))

    for fb, color in zip(fiber_bundles, colors):
        for f in fb:
            if f.shape[0] <= 1:
                continue

            df = f[1:, :] - f[0:-1, :]
            phi = np.arctan2(df[:, 1], df[:, 0])
            theta = np.arccos(df[:, 2] / np.linalg.norm(df, axis=1))
            gphi.extend(phi)
            gtheta.extend(theta)

            if plot == 'scatter':
                phi, theta = _remap_orientation(phi, theta)
                ax.scatter(phi, np.rad2deg(theta), color=color)

    gphi, gtheta = _remap_orientation(np.array(gphi), np.array(gtheta))

    if plot == 'histogram':
        abins = np.linspace(0, 2 * np.pi, 100)
        rbins = np.linspace(0, 90, 50)
        # rbins = np.linspace(0, 0.5 * np.pi, 50)

        #calculate histogram
        hist, _, _ = np.histogram2d(gphi,
                                    np.rad2deg(gtheta),
                                    bins=(abins, rbins))
        A, R = np.meshgrid(abins, rbins)

        pc = ax.pcolormesh(A, R, hist.T, cmap="viridis")
        fig.colorbar(pc)

    ax.set_rmax(90)
    ax.set_rticks(rticks)
    ax.set_rlabel_position(22.5)
    ax.set_yticklabels([])
    ax.grid(True)

    ax.set_title("", va='bottom')
    plt.show()

    return gphi, gtheta
