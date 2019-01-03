import numpy as np

from pyglet.gl import *
from pyglet.window import key
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *

from ..objects import Fiber

WIDTH = 800
HEIGHT = 600
INCREMENT = 5
NUM_POLYGON = 8


class Vertex:
    def __init__(self):
        self.p = np.empty(3, np.float32)
        self.n = np.empty(3, np.float32)


class Window(pyglet.window.Window):
    # Cube 3D start rotation
    xRotation = yRotation = 30
    vertex_lists = []

    def __init__(self, width, height, title=''):
        super(Window, self).__init__(width, height, title)
        glClearColor(0, 0, 0, 1)
        glEnable(GL_DEPTH_TEST)

        glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE)
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)
        glEnable(GL_COLOR_MATERIAL)

        glLightfv(GL_LIGHT0, GL_AMBIENT, (0.2, 0.2, 0.2, 1.0))
        glLightfv(GL_LIGHT0, GL_DIFFUSE, [0.8, 0.8, 0.8, 1.0])
        glLightfv(GL_LIGHT0, GL_SPECULAR, [0.8, 0.8, 0.8, 1.0])

    def on_draw(self):
        self.clear()

        glPushMatrix()

        glRotatef(self.xRotation, 1, 0, 0)
        glRotatef(self.yRotation, 0, 1, 0)

        glColor3f(0.8, 0.80, 0.80)

        for vl in self.vertex_lists:
            vl.draw(pyglet.gl.GL_TRIANGLE_STRIP)

        glPopMatrix()

    def on_resize(self, width, height):
        # set the Viewport
        glViewport(0, 0, width, height)

        # using Projection mode
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()

        aspectRatio = width / height
        gluPerspective(45, aspectRatio, 0.1, 1000000)

        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        glTranslatef(0, 0, -1250)

    def on_text_motion(self, motion):
        if motion == key.UP:
            self.xRotation -= INCREMENT
        elif motion == key.DOWN:
            self.xRotation += INCREMENT
        elif motion == key.LEFT:
            self.yRotation -= INCREMENT
        elif motion == key.RIGHT:
            self.yRotation += INCREMENT


class Vis:
    def __init__(self):
        self.window = Window(WIDTH, HEIGHT, 'fastpli.model.Visualizer')

    def set_data(self, fiber_bundles):
        tube_elm_0 = tuple([Vertex() for _ in range(NUM_POLYGON)])

        # generate circle
        tangent_old = np.array([0, 0, 1], np.float32)
        for k in range(NUM_POLYGON):
            theta = 2 * np.pi * k / NUM_POLYGON
            tube_elm_0[k].n = np.array(
                (np.cos(theta), np.sin(theta), 0), np.float32)

        # convert fb into tubes
        for fb in fiber_bundles:
            for f in fb:
                tube_elm = tube_elm_0[:]
                (points, radii) = f.data
                data_list = np.empty(
                    (len(radii), NUM_POLYGON, 2, 3), np.float32)

                for i in range(points.shape[0]):
                    if i == 0:
                        p0 = points[i, :]
                        p1 = points[i + 1, :]
                        pm = p0
                    elif i == (points.shape[0] - 1):
                        p0 = points[i - 1, :]
                        p1 = points[i, :]
                        pm = p1
                    else:
                        p0 = points[i - 1, :]
                        pm = points[i, :]
                        p1 = points[i + 1, :]

                    if np.array_equal(p0, p1):
                        continue

                    tangent = p1 - p0
                    tangent = tangent / np.linalg.norm(tangent)

                    # rotate old points onto new plane
                    v = np.cross(tangent_old, tangent)
                    s = np.linalg.norm(v) + 1e-9
                    c = np.dot(tangent_old, tangent)
                    rot = np.array(
                        [[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]], np.float32)
                    rot = (np.identity(3, np.float32) + rot) + \
                        np.dot(rot, rot) * (1 - c) / (s * s)

                    for k in range(NUM_POLYGON):
                        p = np.dot(rot, tube_elm[k].n)
                        p = p / np.linalg.norm(p)
                        tube_elm[k].p = pm + p * radii[i]
                        tube_elm[k].n = p

                        data_list[i, k, 0, :] = tube_elm[k].p
                        data_list[i, k, 1, :] = tube_elm[k].n

                c = 0
                vertex_list = pyglet.graphics.vertex_list(
                    (NUM_POLYGON + 1) * (len(radii) - 1) * 2, 'v3f', 'n3f')
                for i in range(len(radii) - 1):
                    for k in range(NUM_POLYGON + 1):
                        vertex_list.vertices[c:c +
                                             3] = data_list[i, k %
                                                            NUM_POLYGON, 0, :]
                        vertex_list.normals[c:c +
                                            3] = data_list[i, k %
                                                           NUM_POLYGON, 1, :]
                        c = c + 3
                        vertex_list.vertices[c:c +
                                             3] = data_list[i +
                                                            1, k %
                                                            NUM_POLYGON, 0, :]
                        vertex_list.normals[c:c +
                                            3] = data_list[i +
                                                           1, k %
                                                           NUM_POLYGON, 1, :]
                        c = c + 3

                self.window.vertex_lists.append(vertex_list)

    def run(self):
        pyglet.app.run()
