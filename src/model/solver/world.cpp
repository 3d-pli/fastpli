#include "world.hpp"

#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <utility>

#include "fiber_class.hpp"
#include "include/aabb.hpp"
#include "include/vemath.hpp"
#include "objects/fiber.hpp"
#include "oct_tree.hpp"

#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num() { return 0; }
inline omp_int_t omp_get_max_threads() { return 1; }
inline omp_int_t omp_get_num_procs() { return 1; }
inline void omp_set_num_threads(int i) {
   (void)i;
   return;
}
#endif

object::fiber::Bundles World::get_fibers() const {
   object::fiber::Bundles fiber_bundles;

   if (fibers_.empty())
      return fiber_bundles;

   fiber_bundles.push_back(std::vector<object::fiber::Fiber>());
   fiber_bundles[0].push_back(fibers_[0]);

   // fiber order does not change, therefore map has the same order
   for (size_t i = 1; i < fibers_.size(); i++) {
      if (map_fb_idx_.at(i - 1).first != map_fb_idx_.at(i).first)
         fiber_bundles.push_back(std::vector<object::fiber::Fiber>());

      fiber_bundles.back().push_back(fibers_[i]);
   }

   return fiber_bundles;
}

void World::set_fibers(
    const std::vector<std::vector<object::fiber::Fiber>> &fiber_bundles) {

   // free memory
   fibers_ = std::vector<geometry::Fiber>();
   map_fb_idx_ = std::map<size_t, std::pair<size_t, size_t>>();

   size_t fb_idx = 0;
   for (auto const &fb : fiber_bundles) {
      size_t f_idx = 0;
      for (auto const &f : fb) {
         map_fb_idx_[fibers_.size()] = std::make_pair(fb_idx, f_idx);
         fibers_.push_back(geometry::Fiber(f, fibers_.size()));
         f_idx++;
      }
      fb_idx++;
   }
}

int World::set_omp_num_threads(int i) {

   if (i > omp_get_num_procs())
      omp_set_num_threads(omp_get_num_procs());
   else
      omp_set_num_threads(i);

   return omp_get_max_threads();
}

bool World::Step() {

   bool solved = true;

   // calculate voi
   aabb::AABB<float, 3> voi{};
   // get first valid voi
   for (auto const &fiber : fibers_) {
      if (!fiber.points().empty()) {
         voi = aabb::AABB<float, 3>(fiber.points()[0]);
         break;
      }
   }
   // combine all vois
   num_obj_ = 0;
   for (auto const &fiber : fibers_) {
      num_obj_ += fiber.ConeSize();
      for (size_t i = 0; i < fiber.size(); i++)
         voi.Unite(aabb::AABB<float, 3>(fiber.points()[i] - fiber.radii()[i],
                                        fiber.points()[i] + fiber.radii()[i],
                                        true));
   }

   // auto smallest_radius = std::numeric_limits<float>::max();
   // auto max_radius = std::numeric_limits<float>::max();
   // for (auto const &fiber : fibers_)
   //    for (auto const &r : fiber.radii)
   //       if (smallest_radius > r)
   //          smallest_radius = r;

   // TODO:
   // if (w_parameter_.obj_mean_length != 0)
   //    if (w_parameter_.obj_mean_length < 3 / 2 * smallest_radius)
   //       w_parameter_.obj_mean_length = 3 / 2 * smallest_radius;

   auto max_obj_size = 3 * w_parameter_.obj_mean_length;
   if (max_obj_size == 0) {
      for (auto const &fiber : fibers_) {
         if (fiber.size() >= 2) {
            for (size_t i = 0; i < fiber.size() - 1; i++) {
               auto delta =
                   vm::length(fiber.points()[i + 1] - fiber.points()[i]);
               delta += std::max(fiber.radii()[i + 1], fiber.radii()[i]);
               if (max_obj_size < delta)
                  max_obj_size = delta;
            }
         }
      }
   }

   // TODO: min_radius sollte gr;-er sein als gleichseitiges dreieck,
   // R=sqrt(3)/3.0*mean_...

   // create OctTree and calculate colliding objects in each leaf
   // TODO: num_threads
   OctTree otree(fibers_, 2 * max_obj_size, col_voi_);
   auto colliding_list = otree.Run();
   num_col_obj_ = colliding_list.size();

   if (!colliding_list.empty()) {
      solved = false;

      // convert set to vector for performance
      std::vector<std::array<size_t, 4>> colliding_vec(colliding_list.begin(),
                                                       colliding_list.end());

      // set speed of colliding objects
#pragma omp parallel for
      for (auto i = 0u; i < colliding_vec.size(); i++) {
         auto elm = colliding_vec[i];

         auto u = fibers_[elm[0]].Cone(elm[1]).PushConesApart(
             fibers_[elm[2]].Cone(elm[3]));
         fibers_[elm[0]].AddSpeed(elm[1], u);
         fibers_[elm[0]].AddSpeed(elm[1] + 1, u);
         fibers_[elm[2]].AddSpeed(elm[3], -u);
         fibers_[elm[2]].AddSpeed(elm[3] + 1, -u);
      }
   }

   // check fiber boundry conditions
#pragma omp parallel for
   for (auto i = 0u; i < fibers_.size(); i++) {
      bool flag_length = fibers_[i].CheckLength(w_parameter_.obj_mean_length);
      bool flag_radius = fibers_[i].CheckRadius(w_parameter_.obj_min_radius);
#pragma omp critical
      solved = solved && flag_length && flag_radius;
   }

   // move colliding objects
   if (!solved) {
#pragma omp parallel for
      for (auto i = 0u; i < fibers_.size(); i++) {
         fibers_[i].Move(w_parameter_.drag);
      }
   }

   return solved;
}

#if _VIS_LIBRARIES
#include "scene.hpp"
void World::DrawScene(float rot_x, float rot_y, float rot_z) {
   if (scene_ == nullptr) {
      char arg0[] = "model.solver";
      char *argv[] = {arg0, nullptr};
      int argc = 1;
      scene_ = std::make_unique<Scene>(argc, argv);
   }
   scene_->SetViewAngle(rot_x, rot_y, rot_z);
   scene_->DrawScene(fibers_);
}
#else
void World::DrawScene(float rot_x, float rot_y, float rot_z) {
   (void)rot_x;
   (void)rot_y;
   (void)rot_z;

   static bool flag = false;

   if (!flag) {
      flag = true;
      std::cout << "No OpenGl detected due build. Deactivating DrawScene()"
                << std::endl;
   }
}

#endif //_VIS_LIBRARIES
