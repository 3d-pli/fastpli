#include <random>
#include <set>
#include <utility>
#include <vector>

// include src libs
#include "aabb.hpp"
#include "fiber_class.hpp"
#include "oct_tree.hpp"
#include "vemath.hpp"

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
   std::vector<std::vector<object::FiberRawData<float>>> get_fibers() const;
   World::WorldParameter get_parameter() const { return w_parameter_; };

   // setter
   void
   set_fibers(std::vector<std::vector<object::FiberRawData<float>>> fibers);
   void set_parameter(World::WorldParameter p) { w_parameter_ = p; };

   // world
   bool Step();
   bool Step(size_t n);
   void RndMoveAll(float value);
   void RndSeed42(bool flag = true) { rnd_seed_42 = flag; };

 private:
   std::vector<object::Fiber> fibers_;
   std::vector<std::pair<size_t, size_t>> fb_idx_2_f_idx_;
   World::WorldParameter w_parameter_;
   bool rnd_seed_42 = false;

   // world functions
   bool CheckRadius();
   bool CheckLength();
};