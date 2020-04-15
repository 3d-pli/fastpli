#include "world.hpp"

#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <utility>

#include "fiber_class.hpp"
#include "include/aabb.hpp"
#include "include/omp.hpp"
#include "include/vemath.hpp"
#include "objects/fiber.hpp"
#include "oct_tree.hpp"
#include "scene.hpp"

object::FiberBundles World::get_fibers() const {
   object::FiberBundles fiber_bundles;

   if (fibers_.empty())
      return fiber_bundles;

   // push_back first element
   size_t i = 0;
   fiber_bundles.push_back(object::FiberBundle());
   fiber_bundles.back().push_back(fibers_[i]);

   // fiber order does not change, therefore map has the same order
   for (i = 1; i < fibers_.size(); i++) {
      // add fiber_bundle if required
      if (map_fb_idx_.at(i - 1).first != map_fb_idx_.at(i).first)
         fiber_bundles.push_back(object::FiberBundle());

      fiber_bundles.back().push_back(fibers_[i]);
   }

   return fiber_bundles;
}

std::vector<std::vector<std::vector<double>>> World::get_fibers_vector() const {
   std::vector<std::vector<std::vector<double>>> fiber_bundles;

   if (fibers_.empty())
      return fiber_bundles;

   // push_back first element
   size_t i = 0;
   fiber_bundles.push_back(std::vector<std::vector<double>>());
   fiber_bundles.back().push_back(fibers_[i].vector());

   // fiber order does not change, therefore map has the same order
   for (i = 1; i < fibers_.size(); i++) {
      // add fiber_bundle if required
      if (map_fb_idx_.at(i - 1).first != map_fb_idx_.at(i).first)
         fiber_bundles.push_back(std::vector<std::vector<double>>());

      fiber_bundles.back().push_back(fibers_[i].vector());
   }

   return fiber_bundles;
}

void World::set_fibers(const object::FiberBundles &fiber_bundles) {

   // free memory
   fibers_ = std::vector<geometry::Fiber>();
   map_fb_idx_ = std::map<size_t, std::pair<size_t, size_t>>();
   max_speed_ = std::numeric_limits<double>::max();

   size_t fb_idx = 0;
   for (auto const &fb : fiber_bundles) {
      size_t f_idx = 0;
      for (auto const &f : fb) {
         map_fb_idx_[fibers_.size()] = std::make_pair(fb_idx, f_idx);
         fibers_.push_back(geometry::Fiber(f, fibers_.size()));
         max_speed_ = std::min(max_speed_, fibers_.back().max_speed());
         f_idx++;
      }
      fb_idx++;
   }

   ResetObjCounter();
}

void World::set_fibers_vector(
    const std::vector<std::vector<std::vector<double>>> &fiber_bundles) {

   // free memory
   fibers_ = std::vector<geometry::Fiber>();
   map_fb_idx_ = std::map<size_t, std::pair<size_t, size_t>>();
   max_speed_ = std::numeric_limits<double>::max();

   size_t fb_idx = 0;
   for (auto const &fb : fiber_bundles) {
      size_t f_idx = 0;
      for (auto const &f : fb) {
         map_fb_idx_[fibers_.size()] = std::make_pair(fb_idx, f_idx);
         fibers_.push_back(geometry::Fiber(object::Fiber(f), fibers_.size()));
         max_speed_ = std::min(max_speed_, fibers_.back().max_speed());
         f_idx++;
      }
      fb_idx++;
   }

   ResetObjCounter();
}

void World::ResetObjCounter() {
   fiber_overlap_ = 0;
   max_level_ = 0;
   num_obj_ = 0;
   num_col_obj_ = 0;
   for (auto &f : fibers_) {
      f.set_max_speed(max_speed_);
      num_obj_ += f.ConeSize();
   }
}

int World::set_omp_num_threads(int i) {

   if (i > omp_get_num_procs())
      omp_set_num_threads(omp_get_num_procs());
   else
      omp_set_num_threads(i);

   return omp_get_max_threads();
}

bool World::ApplyBoundaryConditions(int max_steps) {
   // check fiber boundary conditions

   bool solved = fibers_.empty();
   for (; max_steps >= 0 && !solved; max_steps--) {

#pragma omp parallel for reduction(&& : solved)
      for (auto i = 0u; i < fibers_.size(); i++) {
         bool flag_length =
             fibers_[i].ApplyConeLengthConstrain(w_parameter_.obj_mean_length);
         bool flag_radius =
             fibers_[i].ApplyCurvatureConstrain(w_parameter_.obj_min_radius);
         solved = flag_length && flag_radius;
      }
   }

   num_obj_ = 0;
   for (auto const &fiber : fibers_)
      num_obj_ += fiber.ConeSize();

   return solved;
}

bool World::Step() {

   bool solved = true;

   num_obj_ = 0;
   for (auto const &fiber : fibers_)
      num_obj_ += fiber.ConeSize();

#pragma omp parallel for
   // applying drag before so that velocity is an indicator for colored
   // visualization
   for (auto i = 0u; i < fibers_.size(); i++) {
      fibers_[i].Drag(w_parameter_.drag);
   }

   auto max_obj_size = 0;
   for (auto const &fiber : fibers_) {
      if (fiber.size() >= 2) {
         for (size_t i = 0; i < fiber.size() - 1; i++) {
            auto delta = vm::length(fiber.points()[i + 1] - fiber.points()[i]);
            delta += std::max(fiber.radii()[i + 1], fiber.radii()[i]);
            if (max_obj_size < delta)
               max_obj_size = delta;
         }
      }
   }

   // TODO: min_radius sollte gr;-er sein als gleichseitiges dreieck,
   // R=sqrt(3)/3.0*mean_...

   // create OctTree and calculate colliding objects in each leaf
   // TODO: num_threads
   OctTree otree(fibers_, 1.5 * max_obj_size, col_voi_);
   auto colliding_list = otree.Run();
   // std::cout << "otree.max_level(): " << otree.max_level() << std::endl;
   num_col_obj_ = colliding_list.size();

   if (!colliding_list.empty()) {
      solved = false;

      // convert set to vector for performance
      std::vector<std::array<size_t, 4>> colliding_vec(colliding_list.begin(),
                                                       colliding_list.end());

      // set speed of colliding objects
      fiber_overlap_ = 0;
#pragma omp parallel for reduction(+ : fiber_overlap_)
      for (auto i = 0u; i < colliding_vec.size(); i++) {
         auto elm = colliding_vec[i];

         auto const [f0, f1, f2, f3, dist] =
             fibers_[elm[0]].Cone(elm[1]).PushConesApart(
                 fibers_[elm[2]].Cone(elm[3]));

         // WARNING: not thread safe
         fibers_[elm[0]].AddSpeed(elm[1], f0);
         fibers_[elm[0]].AddSpeed(elm[1] + 1, f1);
         fibers_[elm[2]].AddSpeed(elm[3], f2);
         fibers_[elm[2]].AddSpeed(elm[3] + 1, f3);

         fiber_overlap_ +=
             1 - dist / (std::max(fibers_[elm[0]].radii()[elm[1]],
                                  fibers_[elm[0]].radii()[elm[1] + 1]) +
                         std::max(fibers_[elm[2]].radii()[elm[3]],
                                  fibers_[elm[2]].radii()[elm[3] + 1]));
      }
   }

   // check fiber boundary conditions
#pragma omp parallel for reduction(&& : solved)
   for (auto i = 0u; i < fibers_.size(); i++) {
      bool flag_radius =
          fibers_[i].ApplyCurvatureConstrain(w_parameter_.obj_min_radius);
      solved = solved && flag_radius;
   }

   if (!solved) {
#pragma omp parallel for reduction(&& : solved)
      for (auto i = 0u; i < fibers_.size(); i++) {
         bool flag_length =
             fibers_[i].ApplyConeLengthConstrain(w_parameter_.obj_mean_length);
         solved = solved && flag_length;
      }
   }

   // move colliding objects
   if (!solved) {
#pragma omp parallel for
      for (auto i = 0u; i < fibers_.size(); i++) {
         fibers_[i].Move();
      }
   }

   return solved;
}

void World::DrawScene(double rot_x, double rot_y, double rot_z, bool only_col) {
   if (scene_ == nullptr) {
      char arg0[] = "model.solver";
      char *argv[] = {arg0, nullptr};
      int argc = 1;
      scene_ = std::make_unique<Scene>(argc, argv);
   }
   scene_->SetViewAngle(rot_x, rot_y, rot_z);
   scene_->DrawScene(fibers_, only_col);
}

void World::CloseScene() { scene_->Close(); }
