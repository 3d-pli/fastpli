import fastpli.simulation
import fastpli.analysis
import fastpli.tools
import fastpli.io
import fastpli.objects

import numpy as np
import h5py
import os

import imageio

np.random.seed(42)

FILE_NAME = os.path.abspath(__file__)
FILE_PATH = os.path.dirname(FILE_NAME)
FILE_BASE = os.path.basename(FILE_NAME)

with h5py.File('/tmp/fastpli.example.' + FILE_BASE + '.h5', 'w') as h5f:
    # save script
    h5f['version'] = fastpli.__version__
    with open(os.path.abspath(__file__), 'r') as f:
        h5f['parameter/script'] = f.read()
        h5f['parameter/pip_freeze'] = fastpli.tools.helper.pip_freeze()

    # Setup Simpli for Tissue Generation
    simpli = fastpli.simulation.Simpli()
    simpli.omp_num_threads = 2
    simpli.voxel_size = 1  # in mu meter
    simpli.set_voi([-100, -100, -25], [2350, 550, 25])  # in mu meter

    fiber_bundles = fastpli.io.fiber_bundles.load(
        os.path.join(FILE_PATH, 'fastpli.dat'))
    # fiber_bundles = fastpli.objects.fiber_bundles.Rescale(fiber_bundles, 0.1)
    simpli.fiber_bundles = fiber_bundles

    # (layer_scale, dn, mu, optical-axis)
    # b: background, p: parallel, r: radial
    simpli.fiber_bundles_properties = [[(1, -0.001, 10, 'p')]]

    print('VOI:', simpli.get_voi())
    print('Memory:', str(round(simpli.memory_usage('MB'), 2)) + ' MB')

    # Generate Tissue
    print("Run Generation:")
    label_field, vec_field, tissue_properties = simpli.generate_tissue()

    h5f['tissue/label_field'] = label_field.astype(np.uint16)
    h5f['tissue/vectorfield'] = vec_field
    h5f['tissue/tissue_properties'] = tissue_properties

    # Simulate PLI Measurement
    simpli.filter_rotations = np.deg2rad(np.linspace(0, 180, 18 * 1))
    simpli.light_intensity = 26000  # a.u.
    simpli.interpolate = True
    simpli.wavelength = 525  # in nm
    simpli.resolution = 10  # in mu meter
    simpli.sensor_gain = 3
    simpli.optical_sigma = 0.71  # in voxel size
    # simpli.tilts = np.deg2rad([(0, 0)])
    simpli.tilts = np.deg2rad([(0, 0), (5.5, 0), (5.5, 90), (5.5, 180),
                               (5.5, 270)])

    tilting_stack = [None] * 5
    print("Run Simulation:")
    for t, (theta, phi) in enumerate(simpli.tilts):
        print(round(np.rad2deg(theta), 1), round(np.rad2deg(phi), 1))
        images = simpli.run_simulation(label_field, vec_field,
                                       tissue_properties, theta, phi)

        h5f['simulation/data/' + str(t)] = images

        # apply optic to simulation
        print(images.shape)
        images = simpli.apply_optic(images)
        h5f['simulation/optic/' + str(t)] = images

        # calculate modalities
        epa = simpli.apply_epa(images)
        h5f['analysis/epa/' + str(t) + '/transmittance'] = epa[0]
        h5f['analysis/epa/' + str(t) + '/direction'] = np.rad2deg(epa[1])
        h5f['analysis/epa/' + str(t) + '/retardation'] = epa[2]

        tilting_stack[t] = images

    # save mask for analysis
    mask = np.sum(label_field, 2) > 0
    mask = simpli.apply_optic_resample(1.0 * mask) > 0.1
    h5f['simulation/optic/mask'] = np.uint8(mask)
    mask = None  # keep analysing all pixels

    # print("Run ROFL analysis:")
    rofl_direction, rofl_incl, rofl_t_rel, _ = simpli.apply_rofl(tilting_stack,
                                                                 mask=mask)

    # h5f['analysis/rofl/direction'] = np.rad2deg(rofl_direction)
    # h5f['analysis/rofl/inclination'] = np.rad2deg(rofl_incl)
    # h5f['analysis/rofl/trel'] = rofl_t_rel


    def data2image(data):
        return np.swapaxes(data, 0, 1)
        # return np.swapaxes(np.flip(data, 1), 0, 1)

    imageio.imwrite(os.path.join(FILE_PATH, 'simpli_transmittance.png'),
                    data2image(epa[0]))

    imageio.imwrite(os.path.join(FILE_PATH, 'simpli_direction.png'),
                    data2image(epa[1]))

    imageio.imwrite(os.path.join(FILE_PATH, 'simpli_retardation.png'),
                    data2image(epa[2]))

    imageio.imwrite(os.path.join(FILE_PATH, 'simpli_rofl_direction.png'),
                    data2image(rofl_direction))

    imageio.imwrite(os.path.join(FILE_PATH, 'simpli_rofl_inclination.png'),
                    data2image(rofl_incl))

    imageio.imwrite(os.path.join(FILE_PATH, 'simpli_rofl_trel.png'),
                    data2image(rofl_t_rel))

    imageio.imwrite(
        os.path.join(FILE_PATH, 'simpli_fom.png'),
        data2image(
            fastpli.analysis.images.fom_hsv_black(rofl_direction, rofl_incl)))

    import matplotlib
    import matplotlib.pyplot as plt

    import matplotlib.animation as animation

    fps = 10
    nSeconds = 5
    # snapshots = [np.random.rand(5, 5) for _ in range(nSeconds * fps)]
    snapshots = [data2image(images[:, :, i]) for i in range(images.shape[-1])]

    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure(figsize=(8, 8))

    a = snapshots[0]
    im = plt.imshow(a,
                    interpolation='none',
                    aspect='auto',
                    vmin=np.amin(images),
                    vmax=np.amax(images))
    plt.axis('off')
    plt.axis('equal')
    plt.tight_layout()

    def animate_func(i):
        if i % fps == 0:
            print('.', end='')

        im.set_array(snapshots[i])
        return [im]

    anim = animation.FuncAnimation(
        fig,
        animate_func,
        frames=len(snapshots),
        interval=1000 / fps  # in ms
    )

    anim.save('docs/simpli.gif', writer='imagemagick', fps=fps, dpi=25)

    print('Done!')

    # plt.show()  # Not required, it seems!
