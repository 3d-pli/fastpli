#include <map>
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
   size_t NumObj() const;
   std::vector<std::vector<data::Fiber>> get_fibers() const;
   World::WorldParameter get_parameter() const { return w_parameter_; };

   // setter
   void set_fibers(std::vector<std::vector<data::Fiber>> fibers);
   void set_parameter(World::WorldParameter p) { w_parameter_ = p; };

   // world
   bool Step();
   bool Step(size_t n);

 private:
   std::vector<object::Fiber> fibers_;
   std::map<size_t, std::pair<size_t, size_t>> map_fb_idx_;
   World::WorldParameter w_parameter_;

   // world functions
   bool CheckRadius();
   bool CheckLength();
};