#include "scene.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <GL/gl.h>
#include <GL/glut.h>

#include "include/vemath.hpp"

Scene::Scene(int argc, char **argv) {

   glutInit(&argc, argv);

   glutInitWindowPosition(0, 0);
   glutInitWindowSize(800, 800);
   glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
   glutCreateWindow("fastpli.model.Solver.Visualizer");

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

void Scene::DrawScene(const std::vector<geometry::Fiber> &fibers) {

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
   glTranslatef(-center_.x(), -center_.y(), -center_.z() - 2.5 * distance_);
   glRotated(rotation_.x(), 1, 0, 0);
   glRotated(rotation_.y(), 0, 1, 0);
   glRotated(rotation_.z(), 0, 0, 1);

   DrawCylinders(fibers);

   glutSwapBuffers();
}

void Scene::DrawCylinders(const std::vector<geometry::Fiber> &fibers) {
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
         glTranslatef(points[i].x(), points[i].y(), points[i].z());
         glRotatef(phi, 0.0, 0.0, 1.0);
         glRotatef(theta, 0.0, 1.0, 0.0);
         gluCylinder(quadObj_, radii[i], radii[i + 1], h, 6, 1);
         glPopMatrix();
      }
   }
}

void Scene::AutoVolume(const std::vector<geometry::Fiber> &fibers) {

   vm::Vec3<float> v_min(std::numeric_limits<float>::max());
   vm::Vec3<float> v_max(-std::numeric_limits<float>::max());

   for (const auto &fiber : fibers) {
      for (const auto &pos : fiber.points()) {
         for (int i = 0; i < 3; i++) {
            v_min[i] = std::min(v_min[i], pos[i]);
            v_max[i] = std::max(v_max[i], pos[i]);
         }
      }
   }

   center_new_ = (v_max + v_min) / 2;
   distance_new_ = vm::max(v_max - v_min);

   if ((vm::length(center_ - center_new_) / vm::length(center_)) >
       repos_threshold_)
      center_ = center_new_;
   if (std::abs((distance_ - distance_new_) / distance_) > repos_threshold_)
      distance_ = distance_new_;
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

   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   glViewport(0, 0, w, h);
   gluPerspective(45.0f, ratio, 0.1f, 100000.0f);
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
