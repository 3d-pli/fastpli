#include "scene.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <GL/gl.h>
#include <GL/glut.h>

#include "include/vemath.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

Scene::Scene(int argc, char **argv) {

   glutInit(&argc, argv);

   glutInitWindowPosition(100, 100);
   glutInitWindowSize(800, 800);
   glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
   glutCreateWindow("Vollume Colliding Solver");

   // Lighting set up
   glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   glEnable(GL_COLOR_MATERIAL);

   // Set lighting intensity and color
   GLfloat qaAmbientLight[] = {0.2, 0.2, 0.2, 1.0};
   GLfloat qaDiffuseLight[] = {0.8, 0.8, 0.8, 1.0};
   GLfloat qaSpecularLight[] = {1.0, 1.0, 1.0, 1.0};
   glLightfv(GL_LIGHT0, GL_AMBIENT, qaAmbientLight);
   glLightfv(GL_LIGHT0, GL_DIFFUSE, qaDiffuseLight);
   glLightfv(GL_LIGHT0, GL_SPECULAR, qaSpecularLight);

   glEnable(GL_DEPTH_TEST);

   quadObj_ = gluNewQuadric();
}

Scene::~Scene() {}

void Scene::DrawScene(const std::vector<object::Fiber> &fibers) {

   AutoVolume(fibers);
   CheckWindowSize();

   // set backgound
   glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   glLoadIdentity();

   // set lighning position
   GLfloat light_position[] = {1.0, 1.0, 1.0, 0.0};
   glLightfv(GL_LIGHT0, GL_POSITION, light_position);

   // set axis position
   glTranslatef(offset_.x(), offset_.y(), offset_.z());
   glRotated(rotation_.x(), 1, 0, 0);
   glRotated(rotation_.y(), 0, 1, 0);
   glRotated(rotation_.z(), 0, 0, 1);

   for (const auto &fiber : fibers) {
      if (fiber.size() <= 1)
         continue;

      for (size_t i = 0; i < fiber.size() - 1; i++) {
         auto const &points = fiber.points();
         auto const &radii = fiber.radii();

         auto dp = points[i + 1] - points[i];
         auto h = vm::length(dp);
         auto unit = dp / h;
         auto theta = std::acos(unit.z()) / M_PI * 180;
         auto phi = std::atan2(unit.y(), unit.x()) / M_PI * 180;

         glColor3f(0.8f, 0.8f, 0.8f);
         // TODO: glColor3f(cone.color().x(), cone.color().y(),
         // cone.color().z());

         glPushMatrix();
         glTranslatef(points[i + 1].x(), points[i + 1].y(), points[i + 1].z());
         glRotatef(phi, 0.0, 0.0, 1.0);
         glRotatef(theta, 0.0, 1.0, 0.0);
         gluCylinder(quadObj_, radii[i], radii[i + 1], h, 6, 1);
         glPopMatrix();
      }
   }
   glutSwapBuffers();
}

void Scene::AutoVolume(const vector<object::Fiber> &fibers) {

   vm::Vec3<float> v_min(std::numeric_limits<float>::max());
   vm::Vec3<float> v_max(-std::numeric_limits<float>::max());

   for (const auto &fiber : fibers) {
      for (const auto &pos : fiber.points()) {
         if (pos.x() < v_min.x())
            v_min.x() = pos.x();
         if (pos.y() < v_min.y())
            v_min.y() = pos.y();
         if (pos.z() < v_min.z())
            v_min.z() = pos.z();
         if (pos.x() > v_max.x())
            v_max.x() = pos.x();
         if (pos.y() > v_max.y())
            v_max.y() = pos.y();
         if (pos.z() > v_max.z())
            v_max.z() = pos.z();
      }
   }

   auto volume_dim = v_max - v_min;
   auto center = volume_dim * 0.5f + v_min;

   auto vol_diff = volume_dim - center;
   auto max =
       std::fmax(std::abs(vol_diff.x()),
                 std::fmax(std::abs(vol_diff.y()), std::abs(vol_diff.z())));
   max += v_max.z();

   auto offset = vm::Vec3<float>(0.0f, 0.0f, -max * 1.5f);

   // only change view if nessecary -> smoother video
   if (vm::length(offset_ - offset) > 10) {
      offset_ = offset;
   }
}

void Scene::SetViewAngle(const float x, const float y, const float z) {
   rotation_ = vm::Vec3<float>(x, y, z);
}

void Scene::CheckWindowSize() {

   auto w = glutGet(GLUT_WINDOW_WIDTH);
   auto h = glutGet(GLUT_WINDOW_HEIGHT);

   if (h == 0)
      h = 1;
   float ratio = w * 1.0 / h;

   // Use the Projection Matrix
   glMatrixMode(GL_PROJECTION);

   // Reset Matrix
   glLoadIdentity();

   // Set the viewport to be the entire window
   glViewport(0, 0, w, h);

   // Set the correct perspective.
   gluPerspective(45.0f, ratio, 0.1f, 100000.0f);

   // Get Back to the Modelview
   glMatrixMode(GL_MODELVIEW);
}

void Scene::SavePPM(const char *fname, int start_x, int start_y) {

   auto w = glutGet(GLUT_WINDOW_WIDTH);
   auto h = glutGet(GLUT_WINDOW_HEIGHT);

   FILE *f = fopen(fname, "wb");
   if (!f)
      return;
   std::vector<unsigned char> out(3 * w * h);
   glPixelStorei(GL_PACK_ALIGNMENT, 1); /* byte aligned output */
   glReadPixels(start_x, start_y, w, h, GL_RGB, GL_UNSIGNED_BYTE, &out[0]);
   fprintf(f, "P6\n%d %d\n255\n", w, h);
   for (int y = 0; y < h; y++) { /* flip image bottom-to-top on output */
      fwrite(&out[3 * (h - 1 - y) * w], 1, 3 * w, f);
   }
   fclose(f);
}