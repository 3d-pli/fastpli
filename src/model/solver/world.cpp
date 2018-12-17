#include "world.hpp"

#include <functional>
#include <random>
#include <utility>

#include "fiber_class.hpp"
#include "include/aabb.hpp"
#include "include/vemath.hpp"
#include "objects/fiber.hpp"
#include "oct_tree.hpp"

std::vector<std::vector<data::Fiber>> World::get_fibers() const {
   std::vector<std::vector<data::Fiber>> fiber_bundles;

   auto fb_id_last = fb_idx_2_f_idx_.front().first;
   fiber_bundles.push_back(std::vector<data::Fiber>());
   for (auto const &[fb_id, f_id] : fb_idx_2_f_idx_) {
      if (fb_id != fb_id_last)
         fiber_bundles.push_back(std::vector<data::Fiber>());

      assert(fiber_bundles.size() == fb_id + 1);
      fiber_bundles.back().push_back(fibers_[f_id]);
   }

   return fiber_bundles;
};

void World::set_fibers(std::vector<std::vector<data::Fiber>> fibers) {
   {
      auto tmp = std::vector<object::Fiber>();
      fibers_.swap(tmp);
   }

   fibers_.reserve(fibers.size());
   size_t ind = 0;
   for (size_t i = 0; i < fibers.size(); i++) {
      auto const &fb = fibers[i];
      for (auto const &f : fb) {
         fibers_.push_back(object::Fiber(f, ind));
         fb_idx_2_f_idx_.push_back(std::make_pair(i, ind));
         ind++;
      }
   }
};

bool World::Step() {

   bool solved = true;

   // calculate voi
   aabb::AABB<float, 3> voi{};
   for (auto const &fiber : fibers_) {
      if (!fiber.points().empty()) {
         voi = aabb::AABB<float, 3>(fiber.points()[0]);
         break;
      }
   }

   size_t num_obj = 0;
   for (auto const &fiber : fibers_) {
      num_obj += fiber.ConeSize();
      for (auto const &p : fiber.points()) {
         voi.Unite(p);
      }
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
         num_obj += fiber.ConeSize();
         if (fiber.size() >= 2) {
            for (size_t i = 0; i < fiber.size() - 1; i++) {
               {
                  auto delta =
                      vm::length(fiber.points()[i + 1] - fiber.points()[i]);
                  delta += std::max(fiber.radii()[i + 1], fiber.radii()[i]);
                  if (max_obj_size < delta)
                     max_obj_size = delta;
               }
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

// move colliding objects
#pragma omp parallel for
      for (auto i = 0u; i < fibers_.size(); i++) {
         fibers_[i].Move(w_parameter_.drag);
      }
   }

// check fiber boundry conditions
#pragma omp parallel for
   for (auto i = 0u; i < fibers_.size(); i++) {
      bool flag_length = fibers_[i].CheckLength(w_parameter_.obj_mean_length);
      bool flag_radius = fibers_[i].CheckRadius(w_parameter_.obj_min_radius);
#pragma omp critical
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
