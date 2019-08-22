#include <map>
#include <memory>
#include <set>
#include <utility>
#include <vector>

// include src libs
#include "fiber_class.hpp"
#include "include/aabb.hpp"
#include "include/omp.hpp"
#include "include/vemath.hpp"
#include "objects/fiber.hpp"
#include "oct_tree.hpp"

#if _VIS_LIBRARIES
#include "scene.hpp"
#endif //_VIS_LIBRARIES

class World {
 public:
   struct WorldParameter {
      double drag{0};
      double obj_min_radius{0};
      double obj_mean_length{1};
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
   object::FiberBundles get_fibers() const;
   std::vector<std::vector<std::vector<double>>> get_fibers_vector() const;
   World::WorldParameter get_parameter() const { return w_parameter_; };

   // setter
   int set_omp_num_threads(int i);
   void set_fibers(const object::FiberBundles &fibers);
   void set_fibers_vector(
       const std::vector<std::vector<std::vector<double>>> &fibers);
   void set_parameter(World::WorldParameter p) { w_parameter_ = p; };
   void set_colliding_voi(const aabb::AABB<double, 3> voi) { col_voi_ = voi; };

   // world
   bool Step();
   bool BoundryChecking(int max_steps);
   void DrawScene(double rot_x = 0, double rot_y = 0, double rot_z = 0,
                  bool only_col = false);

 private:
   std::vector<geometry::Fiber> fibers_;
   std::map<size_t, std::pair<size_t, size_t>> map_fb_idx_;
   World::WorldParameter w_parameter_;

   size_t num_obj_{0};
   size_t num_col_obj_{0};

   aabb::AABB<double, 3> col_voi_ = aabb::AABB<double, 3>(vm::Vec3<double>(0));

#if _VIS_LIBRARIES
   std::unique_ptr<Scene> scene_ = nullptr;
#endif //_VIS_LIBRARIES

   // world functions
   bool CheckRadius();
   bool CheckLength();
};
