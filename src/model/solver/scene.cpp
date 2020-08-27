#include "scene.hpp"

#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#elif _WIN32
#include <windows.h>
#else
#include <GL/gl.h>
#include <GL/glut.h>
#endif

#include <GL/freeglut.h>
#include <Python.h>

#include "include/vemath.hpp"

Scene::Scene(int argc, char **argv) {
   argc_ = argc;
   argv_ = new char *[argc + 1];
   for (int i = 0; i <= argc; i++)
      argv_[i] = argv[i];
}

void Scene::Close() {
   if (glut_window_) {
      glutDestroyWindow(glut_window_);
      glut_window_ = 0;
      glutMainLoopEvent();
      glutMainLoopEvent();
   }
}

void Scene::CreateWindow() {

   glutInitWindowPosition(0, 0);
   glutInitWindowSize(800, 800);
   glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
   glutInit(&argc_, argv_);
   quadObj_ = gluNewQuadric();

   glut_window_ = glutCreateWindow("fastpli.model.Solver.Visualizer");

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
}

void Scene::DrawScene(const std::vector<geometry::Fiber> &fibers) {

   if (glut_window_ == 0)
      CreateWindow();

   AutoVolume(fibers);
   CheckWindowSize();

   // set background
   glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   glLoadIdentity();

   // set lightning position
   GLfloat light_position[] = {1.0, 1.0, 1.0, 0.0};
   glLightfv(GL_LIGHT0, GL_POSITION, light_position);

   // set axis position
   glTranslatef(0, 0, -2.5 * distance_);
   glRotated(rotation_.x(), 1, 0, 0);
   glRotated(rotation_.y(), 0, 1, 0);
   glRotated(rotation_.z(), 0, 0, 1);
   glTranslatef(-center_.x(), -center_.y(), -center_.z());

   DrawCylinders(fibers);

   if (axes_)
      DrawAxis();

   glutSwapBuffers();
   glutMainLoopEvent();
}

void Scene::DrawCylinders(const std::vector<geometry::Fiber> &fibers) {
   for (const auto &fiber : fibers) {
      if (fiber.size() <= 1)
         continue;

      for (size_t i = 0; i < fiber.size() - 1; i++) {
         auto const &points = fiber.points();
         auto const &radii = fiber.radii();
         auto const &speed = fiber.speed();

         auto dp = points[i + 1] - points[i];
         auto h = vm::length(dp);
         auto unit = dp / h;
         auto theta = std::acos(unit.z()) / M_PI * 180;
         auto phi = std::atan2(unit.y(), unit.x()) / M_PI * 180;

         glColor3f(0.8f, 0.8f, 0.8f);

         if (speed[i] != vm::Vec3<double>(0) ||
             speed[i + 1] != vm::Vec3<double>(0))
            glColor3f(0.8f, 0.0f, 0.0f);
         else if (only_col_)
            continue;

         glPushMatrix();
         glTranslatef(points[i].x(), points[i].y(), points[i].z());
         glRotatef(phi, 0.0, 0.0, 1.0);
         glRotatef(theta, 0.0, 1.0, 0.0);
         gluCylinder(quadObj_, radii[i], radii[i + 1], h, 6, 1);
         glPopMatrix();
      }
   }
}

void Scene::DrawAxis() {
   glDisable(GL_DEPTH_TEST);

   glLoadIdentity();
   glTranslatef(-center_.x(), -center_.y(), -center_.z() - 2.5 * distance_);
   glRotated(rotation_.x(), 1, 0, 0);
   glRotated(rotation_.y(), 0, 1, 0);
   glRotated(rotation_.z(), 0, 0, 1);

   auto size = distance_ / 5;

   glPushMatrix();
   // x
   glBegin(GL_LINES);
   glColor3f(1.0, 0.0, 0.0);
   glVertex3f(0.0, 0.0f, 0.0f);
   glVertex3f(size, 0.0f, 0.0f);
   glEnd();
   // y
   glBegin(GL_LINES);
   glColor3f(0.0, 1.0, 0.0);
   glVertex3f(0.0, 0.0f, 0.0f);
   glVertex3f(0.0, size, 0.0f);
   glEnd();
   // z
   glBegin(GL_LINES);
   glColor3f(0.0, 0.0, 1.0);
   glVertex3f(0.0, 0.0f, 0.0f);
   glVertex3f(0.0, 0.0f, size);
   glEnd();
   glPopMatrix();

   glEnable(GL_DEPTH_TEST);
}

void Scene::AutoVolume(const std::vector<geometry::Fiber> &fibers) {

   static bool first_time = true;

   vm::Vec3<double> v_min(std::numeric_limits<float>::max());
   vm::Vec3<double> v_max(-std::numeric_limits<float>::max());

   for (const auto &fiber : fibers) {
      for (const auto &pos : fiber.points()) {
         for (int i = 0; i < 3; i++) {
            v_min[i] = std::min(v_min[i], pos[i]);
            v_max[i] = std::max(v_max[i], pos[i]);
         }
      }
   }

   auto center_new_ = vm::cast<float>(v_max + v_min) / 2;
   auto distance_new_ = vm::max(v_max - v_min);
   if ((vm::length(center_ - center_new_) / vm::length(center_)) >
       repos_threshold_)
      center_ = center_new_;
   if (std::abs((distance_ - distance_new_) / distance_) > repos_threshold_)
      distance_ = distance_new_;

   if (first_time) {
      first_time = false;
      center_ = center_new_;
      distance_ = distance_new_;
   }

   if (distance_user_ > 0)
      distance_ = distance_user_;
   if (center_user_)
      center_ = center_user_value_;
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
