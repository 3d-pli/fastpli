#ifndef SRC_MODEL_SOLVER_OCT_TREE_HPP_
#define SRC_MODEL_SOLVER_OCT_TREE_HPP_

#include <array>
#include <set>
#include <tuple>
#include <vector>

#include "fiber_class.hpp"
#include "fiber_segment.hpp"
#include "include/aabb.hpp"
#include "include/vemath.hpp"

class OctTree {
 public:
   OctTree(const std::vector<geometry::Fiber> &fibers,
           const double min_cube_size, const aabb::AABB<double, 3> col_voi);
   std::set<std::array<size_t, 4>> Run();

   int max_level() { return max_level_; }

 private:
   std::tuple<std::vector<std::vector<size_t>>, int>
   GenerateLeafs(const std::vector<size_t> &ids,
                 const aabb::AABB<double, 3> &cube, int level);
   std::vector<std::array<size_t, 4>>
   TestCollision(const std::vector<size_t> &fs_ids);

   std::vector<geometry::FiberSegment> fiber_segments_;
   int max_level_ = 0;
   double min_cube_size_ = 0;
   aabb::AABB<double, 3> voi_cube_{};
   const size_t kMaxParticle_ = 10;
   const int kMaxThreadLevel_ = 1;
};

#endif // SRC_MODEL_SOLVER_OCT_TREE_HPP_
