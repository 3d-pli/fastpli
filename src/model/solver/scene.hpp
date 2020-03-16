#ifndef SCENE_HPP_
#define SCENE_HPP_

#include <vector>

#include "fiber_class.hpp"
#include "include/vemath.hpp"

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#else
#ifdef _WIN32
  #include <windows.h>
#endif
#include <GL/gl.h>
#include <GL/glut.h>
#endif

class Scene {
 public:
   Scene(int argc, char **argv);
   ~Scene() = default;

   void SetViewAngle(const float x, const float y, const float z);
   void DrawScene(const std::vector<geometry::Fiber> &fibers,
                  const bool only_col = false);
   void SavePPM(const char *fname, int start_x = 0, int start_y = 0);

 private:
   void AutoVolume(const std::vector<geometry::Fiber> &fibers);
   void DrawCylinders(const std::vector<geometry::Fiber> &fibers);
   void CheckWindowSize();

   GLUquadricObj *quadObj_;
   bool only_col_ = false;
   vm::Vec3<float> rotation_ = 0;
   vm::Vec3<float> center_ = 0;
   vm::Vec3<float> center_new_ = center_;
   float distance_ = 0;
   float distance_new_ = 0;
   const float repos_threshold_ = 0.25;
};

#endif // SCENECLASS_HPP_
