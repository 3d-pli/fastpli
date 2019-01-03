#include <map>
#include <memory>
#include <random>
#include <set>
#include <utility>
#include <vector>

// include src libs
#include "fiber_class.hpp"
#include "include/aabb.hpp"
#include "include/vemath.hpp"
#include "objects/fiber.hpp"
#include "oct_tree.hpp"
#include "scene.hpp"

class World {
 public:
   struct WorldParameter {
      float drag{0};
      float obj_min_radius{0};
      float obj_mean_length{1};
   };

   // defaults
   World() = default;
   World(World &&) = default;
   World(const World &) = default;
   World &operator=(World &&) = default;
   World &operator=(const World &) = default;
   ~World() = default;

   // getter
   size_t NumObj() const { return num_obj_; };
   size_t NumColObj() const { return num_col_obj_; };
   std::vector<std::vector<data::Fiber>> get_fibers() const;
   World::WorldParameter get_parameter() const { return w_parameter_; };

   // setter
   int set_omp_num_threads(int i);
   void set_fibers(const std::vector<std::vector<data::Fiber>> &fibers);
   void set_parameter(World::WorldParameter p) { w_parameter_ = p; };
   void set_colliding_voi(const aabb::AABB<float, 3> voi) { col_voi_ = voi; };

   // world
   bool Step();
   void Visualization(bool flag) { do_vis = flag; };

 private:
   std::vector<object::Fiber> fibers_;
   std::map<size_t, std::pair<size_t, size_t>> map_fb_idx_;
   World::WorldParameter w_parameter_;

   size_t num_obj_{0};
   size_t num_col_obj_{0};

   bool do_vis = false;
   std::unique_ptr<Scene> scene_ = nullptr;

   aabb::AABB<float, 3> col_voi_ = aabb::AABB<float, 3>(vm::Vec3<float>(0));

   // world functions
   bool CheckRadius();
   bool CheckLength();
};
