#include "world.hpp"

#include <functional>
#include <map>
#include <random>
#include <utility>

#include "fiber_class.hpp"
#include "include/aabb.hpp"
#include "include/vemath.hpp"
#include "objects/fiber.hpp"
#include "oct_tree.hpp"

const std::vector<std::vector<data::Fiber>> World::get_fibers() const {
   std::vector<std::vector<data::Fiber>> fiber_bundles;

   if (fibers_.empty())
      return fiber_bundles;

   fiber_bundles.push_back(std::vector<data::Fiber>());
   fiber_bundles[0].push_back(fibers_[0]);

   // since fiber order did not change
   for (size_t i = 1; i < fibers_.size(); i++) {
      if (map_fb_idx_.at(i - 1).first != map_fb_idx_.at(i).first)
         fiber_bundles.push_back(std::vector<data::Fiber>());

      fiber_bundles.back().push_back(fibers_[i]);
   }

   return fiber_bundles;
}

void World::set_fibers(
    const std::vector<std::vector<data::Fiber>> &fiber_bundles) {

   // TODO: check reset everything
   {
      auto tmp = std::vector<object::Fiber>();
      fibers_.swap(tmp);
   }
   map_fb_idx_.clear();

   size_t fb_idx = 0;
   for (auto const &fb : fiber_bundles) {
      size_t f_idx = 0;
      for (auto const &f : fb) {
         map_fb_idx_[fibers_.size()] = std::make_pair(fb_idx, f_idx);
         fibers_.push_back(object::Fiber(f, fibers_.size()));
         f_idx++;
      }
      fb_idx++;
   }
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
   for (auto const &fiber : fibers_) {
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

   // create OctTree and calculate colliding objects
   // TODO: num_threads
   OctTree otree(fibers_, 2 * max_obj_size);

   // std::cout << "Step(): starting otree: " << 2 * max_obj_size << "\n";
   auto colliding_list = otree.Run();

   // std::cout << "Step(): finished otree: " << colliding_list.size() << "\n";

   if (!colliding_list.empty()) {
      solved = false;

      // convert set to vector for performance
      std::vector<std::array<size_t, 4>> colliding_vec(colliding_list.begin(),
                                                       colliding_list.end());

      // set speed of colliding objects
      // #pragma omp parallel for
      for (auto i = 0u; i < colliding_vec.size(); i++) {
         auto elm = colliding_vec[i];

         auto u = fibers_[elm[0]].Cone(elm[1]).PushConesApart(
             fibers_[elm[2]].Cone(elm[3]));
         fibers_[elm[0]].AddSpeed(elm[1], u);
         fibers_[elm[0]].AddSpeed(elm[1] + 1, u);
         fibers_[elm[2]].AddSpeed(elm[3], -u);
         fibers_[elm[2]].AddSpeed(elm[3] + 1, -u);
      }

      // move colliding objects
      // #pragma omp parallel for
      for (auto i = 0u; i < fibers_.size(); i++) {
         fibers_[i].Move(w_parameter_.drag);
      }
   }

   // check fiber boundry conditions
   // #pragma omp parallel for
   for (auto i = 0u; i < fibers_.size(); i++) {
      bool flag_length = fibers_[i].CheckLength(w_parameter_.obj_mean_length);
      bool flag_radius = fibers_[i].CheckRadius(w_parameter_.obj_min_radius);
      // #pragma omp critical
      solved = solved && flag_length && flag_radius;
      // std::cout << "i:" << i << "\t" << solved << "\t" << flag_length << "\t"
      //           << flag_radius << "\t" << num_obj << "\n";
   }

   return solved;
}

bool World::Step(size_t n) {
   bool flag = false;
   for (size_t i = 0; i < n; i++)
      if ((flag = Step()))
         return true;
   return flag;
}

size_t World::NumObj() const {
   size_t n{0};
   for (auto const &fiber : fibers_)
      n += fiber.ConeSize();
   return n;
}
