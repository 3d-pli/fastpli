import numpy as np
from fastpli.simulation import Simpli

'''
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 1)

class IndexTracker(object):
    def __init__(self, ax, X):
        self.ax = ax
        ax.set_title('use scroll wheel to navigate images')

        self.X = X
        rows, cols, self.slices = X.shape
        self.ind = self.slices//2

        self.im = ax.imshow(self.X[:, :, self.ind])
        self.update()

    def onscroll(self, event):
        # print("%s %s" % (event.button, event.step))
        if event.button == 'up':
            self.ind = (self.ind + 1) % self.slices
        else:
            self.ind = (self.ind - 1) % self.slices
        self.update()

    def update(self):
        self.im.set_data(self.X[:, :, self.ind])
        ax.set_ylabel('slice %s' % self.ind)
        self.im.axes.figure.canvas.draw()
'''

simpli = Simpli()
# PliGeneration ###
simpli.pixel_size = 1
simpli.dim_global = [100, 100, 100]
simpli.dim_local = [100, 100, 100]
simpli.dim_offset_local = [0, 0, 0]

simpli.ReadFiberFile('example/cube.h5')
simpli.SetFiberProperties([[(0.333, 0.004, 10, 'p'), (
    0.666, -0.004, 5, 'b'), (1.0, 0.004, 1, 'r')]])

# manipulation of fibers
simpli.RotateVolumeAroundPoint(np.deg2rad(
    20), np.deg2rad(-10), np.deg2rad(5), [10, -5, 7.5])
simpli.TranslateVolume([25, -15, 50])

label_field, _, tissue_properties = simpli.GenerateTissue(only_label=True)
label_field, vec_field, tissue_properties = simpli.GenerateTissue()

# TODO:
# label_field_vis = label_field.copy()
# label_field_vis[label_field > 0] += 3

# PliSimulation ###
simpli.setup.filter_rotations = np.deg2rad([0, 30, 60, 90, 120, 150])
simpli.setup.light_intensity = 26000
simpli.setup.resolution = 1
simpli.setup.untilt_sensor = True
simpli.setup.wavelength = 525

# simpli.InitSimulation(label_field, vec_field, tissue_properties)

print("run_simulation: 0")
image = simpli.run_simulation(label_field, vec_field, tissue_properties, 0, 0)

print("run_simulation: 1")
image = simpli.run_simulation(
    label_field,
    vec_field,
    tissue_properties,
    np.deg2rad(5.5),
    np.deg2rad(0))

print("run_simulation: 2")
image = simpli.run_simulation(
    label_field,
    vec_field,
    tissue_properties,
    np.deg2rad(5.5),
    np.deg2rad(90))

print("run_simulation: 3")
image = simpli.run_simulation(
    label_field,
    vec_field,
    tissue_properties,
    np.deg2rad(5.5),
    np.deg2rad(180))

print("run_simulation: 4")
image = simpli.run_simulation(
    label_field,
    vec_field,
    tissue_properties,
    np.deg2rad(5.5),
    np.deg2rad(270))

'''
tracker = IndexTracker(ax,image)
fig.canvas.mpl_connect('scroll_event', tracker.onscroll)
plt.show()
'''
