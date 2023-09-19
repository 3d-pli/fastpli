#ifndef SRC_MODEL_SOLVER_WORLD_HPP_
#define SRC_MODEL_SOLVER_WORLD_HPP_

#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

// include src libs
#include "fiber_class.hpp"
#include "include/aabb.hpp"
#include "include/omp.hpp"
#include "include/vemath.hpp"
#include "objects/fiber.hpp"
#include "oct_tree.hpp"

#if _VIS_LIBRARIES
#include "scene.hpp"
#endif //_VIS_LIBRARIES

class World {
 public:
   struct WorldParameter {
      double drag{0};
      double obj_min_radius{0};
      double obj_mean_length{0};
   };

   // defaults
   World() = default;
   World(World &&) = default;
   World &operator=(World &&) = default;
   ~World() = default;

   // getter
   int64_t num_obj() const { return num_obj_; }
   int64_t num_col_obj() const { return num_col_obj_; }
   double overlap() const { return fiber_overlap_; }
   double max_speed() const { return max_speed_; }
   object::FiberBundles get_fibers() const;
   std::vector<std::vector<std::vector<double>>> get_fibers_vector() const;
   World::WorldParameter get_parameter() const { return w_parameter_; }

   // setter
   int set_omp_num_threads(int i);
   void set_fibers(const object::FiberBundles &fibers);
   void set_fibers_vector(
       const std::vector<std::vector<std::vector<double>>> &fibers);
   void set_parameter(World::WorldParameter p) { w_parameter_ = p; }
   void set_colliding_voi(const aabb::AABB<double, 3> voi) { col_voi_ = voi; }

   // world
   bool Step();
   bool ApplyBoundaryConditions(int max_steps);

   // scene
   void InitScene();
   void DrawScene();
   void ResetView() {
#if _VIS_LIBRARIES
      InitScene();
      if (scene_)
         scene_->ResetView();
#endif
   }
   void ToggleAxis(bool flag) {
#if _VIS_LIBRARIES
      InitScene();
      if (scene_)
         scene_->ToggleAxis(flag);
#endif
      (void)flag;
   }

   void ToggleCollisionView(bool flag) {
#if _VIS_LIBRARIES
      InitScene();
      if (scene_)
         scene_->ToggleCollisionView(flag);
#endif
      (void)flag;
   }
   void SetViewCenter(const float x, const float y, const float z) {
#if _VIS_LIBRARIES
      InitScene();
      if (scene_)
         scene_->SetViewCenter(x, y, z);
#endif
      (void)x;
      (void)y;
      (void)z;
   }
   void SetViewAngles(const float x, const float y, const float z) {
#if _VIS_LIBRARIES
      InitScene();
      if (scene_)
         scene_->SetViewAngles(x, y, z);
#endif
      (void)x;
      (void)y;
      (void)z;
   }
   void SetViewDistance(const float d) {
#if _VIS_LIBRARIES
      InitScene();
      if (scene_)
         scene_->SetViewDistance(d);
#endif
      (void)d;
   }
   void CloseScene() {
#if _VIS_LIBRARIES
      InitScene();
      if (scene_ != nullptr)
         scene_->Close();
#endif
   }
   void SavePPM(std::string file) {
#if _VIS_LIBRARIES
      InitScene();
      if (scene_)
         scene_->SavePPM(file.c_str(), 0, 0);
#else
      (void)file;
#endif
   }

   void SaveSTL(const char *fname);

 private:
   std::vector<geometry::Fiber> fibers_;
   std::map<size_t, std::pair<size_t, size_t>> map_fb_idx_;
   aabb::AABB<double, 3> col_voi_{};
   World::WorldParameter w_parameter_;

   double max_speed_{0};
   double fiber_overlap_{-1};

   int64_t max_level_{0};
   int64_t num_obj_{0};
   int64_t num_col_obj_{0};

   bool draw_axis_{false};
#if _VIS_LIBRARIES
   std::unique_ptr<Scene> scene_ = nullptr;
#endif //_VIS_LIBRARIES

   // world functions
   bool ApplyCurvatureConstrain();
   bool ApplyFiberSegmentLengthConstrain();
   void ResetObjCounter();
};

#endif // SRC_MODEL_SOLVER_WORLD_HPP_
