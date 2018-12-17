#ifndef OCT_TREE_HPP_
#define OCT_TREE_HPP_

#include <array>
#include <set>
#include <vector>

#include "cone_class.hpp"
#include "fiber_class.hpp"
#include "include/aabb.hpp"
#include "include/vemath.hpp"

class OctTree {
 public:
   OctTree(const std::vector<object::Fiber> fibers, const float min_cube_size);
   std::set<std::array<size_t, 4>> Run();

 private:
   std::vector<std::vector<size_t>>
   GenerateLeafs(const std::vector<size_t> &ids,
                 const aabb::AABB<float, 3> &cube, int level);
   std::set<std::array<size_t, 4>>
   TestCollision(const std::vector<size_t> &cone_ids);

   std::vector<object::Cone> cones_;
   // TODO: int num_threads_ = 0;
   float min_cube_size_ = 0;
   aabb::AABB<float, 3> main_cube_;
   static constexpr long long max_particle_ = 10;
};

#endif // OCT_TREE_HPP_