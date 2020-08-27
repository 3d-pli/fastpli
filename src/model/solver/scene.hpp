#ifndef SCENE_HPP_
#define SCENE_HPP_

#include <vector>

#include "fiber_class.hpp"
#include "include/vemath.hpp"

#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#elif _WIN32
#include <windows.h>
#else
#include <GL/gl.h>
#include <GL/glut.h>
#endif

class Scene {
 public:
   Scene(int argc, char **argv);
   ~Scene() { delete[] argv_; };

   void ResetView() {
      center_user_ = false;
      distance_user_ = 0;
      rotation_ = {30, 30, 0};
      axes_ = false;
      only_col_ = false;
   };
   void SetViewAngles(const float x, const float y, const float z) {
      rotation_ = {x, y, z};
   };
   void SetViewCenter(const float x, const float y, const float z) {
      center_user_ = true;
      center_user_value_ = {x, y, z};
   };
   void SetViewDistance(const float d) { distance_user_ = d; };
   void DrawScene(const std::vector<geometry::Fiber> &fibers);
   void ToggleAxis(bool flag) { axes_ = flag; };
   void ToggleCollisionView(bool flag) { only_col_ = flag; };
   void SetAxis(bool flag) { axes_ = flag; };
   void SavePPM(const char *fname, int start_x = 0, int start_y = 0);
   void Close();

 private:
   void CreateWindow();
   void AutoVolume(const std::vector<geometry::Fiber> &fibers);
   void DrawCylinders(const std::vector<geometry::Fiber> &fibers);
   void DrawAxis();
   void CheckWindowSize();

   int glut_window_ = 0;

   int argc_ = 0;
   char **argv_;

   GLUquadricObj *quadObj_ = nullptr;

   vm::Vec3<float> rotation_ = {30, 30, 0};
   vm::Vec3<float> center_ = 0;
   vm::Vec3<float> center_user_value_ = 0;

   float distance_ = 0;
   float distance_user_ = 0;
   const float repos_threshold_ = 0.25;

   bool only_col_ = false;
   bool center_user_ = false;
   bool axes_ = false;
};

#endif // SCENECLASS_HPP_
