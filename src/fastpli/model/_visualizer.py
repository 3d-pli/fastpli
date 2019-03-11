import numpy as np

# package PyOpenGL:
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *

NUM_POLYGON = 6
DISTANCE_FACTOR = 2.5


class Visualizer:
    """
    None interactive visualization of FiberBundles.

    Note:
        The methode fastpli.model.Solver.DrawScene() is much faster.
        Later this method will be an interactive visualizer
    """

    fbs = []
    rot_x = rot_y = 30
    distance = 0
    distance_new = distance
    center = np.empty((3), np.float32)
    center_new = center

    def __init__(self, width=800, height=600, title='fastpli.model.Visualizer'):
        self._glut_init(width, height, title)

    def _glut_init(self, width, height, title):
        glutInit()
        glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH)
        glutInitWindowSize(width, height)
        glutInitWindowPosition(0, 0)
        glutCreateWindow(title)

        # setting scene
        glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE)
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)
        glEnable(GL_COLOR_MATERIAL)

        glLightfv(GL_LIGHT0, GL_AMBIENT, (0.2, 0.2, 0.2, 1.0))
        glLightfv(GL_LIGHT0, GL_DIFFUSE, [0.8, 0.8, 0.8, 1.0])
        glLightfv(GL_LIGHT0, GL_SPECULAR, [0.8, 0.8, 0.8, 1.0])

        glEnable(GL_DEPTH_TEST)

    def _draw_cylinders(self):
        quadObj = gluNewQuadric()
        glColor3f(0.8, 0.8, 0.8)
        for fb in self.fbs:
            for f in fb:
                (points, radii) = f.data
                for i in range(len(radii) - 1):
                    dp = points[i + 1, :] - points[i, :]
                    height = np.linalg.norm(dp)
                    unit = dp / height
                    theta = np.arccos(unit[2]) / np.pi * 180
                    phi = np.arctan2(unit[1], unit[0]) / np.pi * 180

                    glPushMatrix()
                    glTranslatef(points[i, 0], points[i, 1], points[i, 2])
                    glRotatef(phi, 0.0, 0.0, 1.0)
                    glRotatef(theta, 0.0, 1.0, 0.0)
                    gluCylinder(quadObj, radii[i], radii[i + 1], height,
                                NUM_POLYGON, 1)
                    glPopMatrix()

    def _auto_volume(self):
        max_vol = np.array((-float('inf'), -float('inf'), -float('inf')))
        min_vol = np.array((float('inf'), float('inf'), float('inf')))

        for fb in self.fbs:
            for f in fb:
                fiber_min = f.points.min(axis=0)
                fiber_max = f.points.max(axis=0)
                min_vol = np.array([min_vol, fiber_min]).min(axis=0)
                max_vol = np.array([max_vol, fiber_max]).max(axis=0)

        self.center_new = (min_vol + max_vol) / 2
        self.distance_new = max(max_vol - min_vol)

    def _resize(self, width, height):
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        glViewport(0, 0, width, height)
        gluPerspective(45, width / height, 0.1, 1000000)
        glMatrixMode(GL_MODELVIEW)

    def set_fbs(self, fiber_bundles):
        self.fbs = fiber_bundles
        self._auto_volume()

        if sum(abs(self.center - self.center_new) / self.center) > 0.1:
            self.center = self.center_new
        if abs((self.distance - self.distance_new) / self.distance) > 0.1:
            self.distance = self.distance_new

    def draw(self):
        self._resize(glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT))

        glClearColor(0, 0, 0, 1)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glLightfv(GL_LIGHT0, GL_POSITION, [1, 1, 1, 0])

        glLoadIdentity()

        glTranslatef(-self.center[0], -self.center[1], -self.center[2])
        glTranslatef(0, 0, -DISTANCE_FACTOR * self.distance)

        glRotatef(self.rot_x, 1, 0, 0)
        glRotatef(self.rot_y, 0, 1, 0)

        self._draw_cylinders()

        glutSwapBuffers()

    def set_rot(self, rotx, roty):
        self.rot_x += rotx
        self.rot_y += roty
