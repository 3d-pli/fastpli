#include "oct_tree.hpp"

#include <algorithm>
#include <array>
#include <tuple>
#include <vector>

#include "cone_class.hpp"
#include "include/aabb.hpp"
#include "include/vemath.hpp"

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

// aabb helper
std::array<aabb::AABB<float, 3>, 8> SplitInto8Cubes(aabb::AABB<float, 3> cube) {

   auto x_max_ = cube.max.x();
   auto y_max_ = cube.max.y();
   auto z_max_ = cube.max.z();
   auto x_min_ = cube.min.x();
   auto y_min_ = cube.min.y();
   auto z_min_ = cube.min.z();
   std::array<aabb::AABB<float, 3>, 8> result;

   auto dx = (x_max_ - x_min_) * 0.5f;
   auto dy = (y_max_ - y_min_) * 0.5f;
   auto dz = (z_max_ - z_min_) * 0.5f;

   vm::Vec3<float> a, b;

   a = {x_min_, y_min_, z_min_};
   b = {x_min_ + dx, y_min_ + dy, z_min_ + dz};
   result[0] = aabb::AABB<float, 3>(a, b);

   a = {x_min_ + dx, y_min_, z_min_};
   b = {x_max_, y_min_ + dy, z_min_ + dz};
   result[1] = aabb::AABB<float, 3>(a, b);

   a = {x_min_, y_min_ + dy, z_min_};
   b = {x_min_ + dx, y_max_, z_min_ + dz};
   result[2] = aabb::AABB<float, 3>(a, b);

   a = {x_min_ + dx, y_min_ + dy, z_min_};
   b = {x_max_, y_max_, z_min_ + dz};
   result[3] = aabb::AABB<float, 3>(a, b);

   a = {x_min_, y_min_, z_min_ + dz};
   b = {x_min_ + dx, y_min_ + dy, z_max_};
   result[4] = aabb::AABB<float, 3>(a, b);

   a = {x_min_ + dx, y_min_, z_min_ + dz};
   b = {x_max_, y_min_ + dy, z_max_};
   result[5] = aabb::AABB<float, 3>(a, b);

   a = {x_min_, y_min_ + dy, z_min_ + dz};
   b = {x_min_ + dx, y_max_, z_max_};
   result[6] = aabb::AABB<float, 3>(a, b);

   a = {x_min_ + dx, y_min_ + dy, z_min_ + dz};
   b = {x_max_, y_max_, z_max_};
   result[7] = aabb::AABB<float, 3>(a, b);

   return result;
}

OctTree::OctTree(const std::vector<geometry::Fiber> &fibers,
                 const float min_cube_size,
                 const aabb::AABB<float, 3> col_voi) {

   min_cube_size_ = min_cube_size;

   // create vector of cones
   for (auto const &fiber : fibers) {
      auto const fiber_cones = fiber.Cones();
      cones_.insert(cones_.end(), std::make_move_iterator(fiber_cones.begin()),
                    std::make_move_iterator(fiber_cones.end()));
   }

   // calculate bounding box
   if (col_voi.min == col_voi.max) {
      // check the whole volume
      for (auto const &fiber : fibers) {
         if (!fiber.points().empty()) {
            main_cube_ = aabb::AABB<float, 3>(fiber.points()[0]);
            break;
         }
      }

      for (auto const &fiber : fibers) {
         for (auto const &p : fiber.points()) {
            main_cube_.Unite(p);
         }
      }
   } else {
      main_cube_ = col_voi;
   }
}

std::set<std::array<size_t, 4>> OctTree::Run() {
   std::set<std::array<size_t, 4>> results;
   std::vector<size_t> ids, all_ids;
   all_ids.resize(cones_.size());
   std::iota(all_ids.begin(), all_ids.end(), 0);

   for (auto id : all_ids) {
      // TODO: test impact of aabb::Overlap
      if (aabb::Overlap(main_cube_, cones_[id].aabb()))
         ids.push_back(id);
   }

   auto const leafs = GenerateLeafs(ids, main_cube_, 0);

#pragma omp parallel for
   for (size_t i = 0; i < leafs.size(); i++) {
      auto const result = TestCollision(leafs[i]);

#pragma omp critical
      results.insert(std::make_move_iterator(result.begin()),
                     std::make_move_iterator(result.end()));
   }

   return results;
}

std::vector<std::vector<size_t>>
OctTree::GenerateLeafs(const std::vector<size_t> &ids,
                       const aabb::AABB<float, 3> &cube, int level) {

   std::vector<std::vector<size_t>> tree_ids;

   if (ids.size() < max_particle_ ||
       (cube.max.x() - cube.min.x()) < 2 * min_cube_size_) {
      tree_ids.push_back(ids);
   } else {

      // TODO: if (level < ceil(log(num_threads_) / log(8))) {
      if (level <= 1) {

         auto sub_cubes = SplitInto8Cubes(cube);

#pragma omp parallel for schedule(static)
         for (auto i = 0; i < 8; i++) {
            std::vector<size_t> sub_ids;
            for (auto id : ids) {
               if (aabb::Overlap(sub_cubes[i], cones_[id].aabb()))
                  sub_ids.push_back(id);
            }

            auto sub_tree = GenerateLeafs(sub_ids, sub_cubes[i], level + 1);

#pragma omp critical
            {
               tree_ids.insert(tree_ids.end(),
                               std::make_move_iterator(sub_tree.begin()),
                               std::make_move_iterator(sub_tree.end()));
            }
         }
      } else {

         auto sub_cubes = SplitInto8Cubes(cube);

         for (auto i = 0; i < 8; i++) {
            std::vector<size_t> sub_ids;
            for (auto id : ids) {
               if (aabb::Overlap(sub_cubes[i], cones_[id].aabb()))
                  sub_ids.push_back(id);
            }

            auto sub_tree = GenerateLeafs(sub_ids, sub_cubes[i], level + 1);
            tree_ids.insert(tree_ids.end(),
                            std::make_move_iterator(sub_tree.begin()),
                            std::make_move_iterator(sub_tree.end()));
         }
      }
   }

   return tree_ids;
}

std::set<std::array<size_t, 4>>
OctTree::TestCollision(const std::vector<size_t> &cone_ids) {
   std::set<std::array<size_t, 4>> result;

   if (cone_ids.empty())
      return result;

   for (auto i = 0u; i < cone_ids.size() - 1; i++) {
      for (auto j = i + 1; j < cone_ids.size(); j++) {
         auto const cone1 = cones_[cone_ids[i]];
         auto const cone2 = cones_[cone_ids[j]];

         if (cone1.fiber_idx == cone2.fiber_idx) {
            // TODO: optimize(static_cast<ll>, all var ll, ...)
            auto const delta = cone1.cone_idx >= cone2.cone_idx
                                   ? cone1.cone_idx - cone2.cone_idx
                                   : cone2.cone_idx - cone1.cone_idx;

            // direct neighbor
            if (delta <= 1)
               continue;

            auto const mean_seg_length =
                vm::length((cone1.p1 - cone1.p0) + (cone2.p1 - cone2.p0)) * 0.5;

            // check for linked chain intersection
            auto const sum_r = cone1.r + cone2.r;

            // TODO: assumption radii did not change much
            if (2.0 / 3.0 * mean_seg_length * delta < sum_r * 2) // 2 for safety
               continue;
         }
         if (cone1.CollideWith(cone2)) {
            result.insert({cone1.fiber_idx, cone1.cone_idx, cone2.fiber_idx,
                           cone2.cone_idx});
         }
      }
   }
   return result;
}
