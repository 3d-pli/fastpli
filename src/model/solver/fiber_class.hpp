#ifndef FIBER_CLASS_HPP_
#define FIBER_CLASS_HPP_

#include <vector>

#include "cone_class.hpp"
#include "include/vemath.hpp"
#include "objects/fiber.hpp"

namespace geometry {
class Fiber : public object::Fiber {
 public:
   // defaults
   Fiber(Fiber &&) = default;
   Fiber(const Fiber &) = default;
   Fiber &operator=(Fiber &&) = default;
   Fiber &operator=(const Fiber &) = default;
   ~Fiber() = default;

   // constructors
   Fiber(const object::Fiber &fiber, const size_t fiber_idx);

   // getter
   const size_t &fiber_idx() const { return fiber_idx_; }
   const std::vector<vm::Vec3<double>> &speed() const { return speed_; }

   // cones
   size_t ConeSize() const;
   object::Cone Cone(size_t i) const;
   std::vector<object::Cone> Cones() const;

   // manipulators
   void Move();
   void Drag(const double drag = 1);
   bool CheckRadius(const double obj_min_radius);
   bool CheckLength(const double obj_mean_length);

   void Split(size_t idx);
   void Combine(size_t idx);

   void AddSpeed(size_t idx, const vm::Vec3<double> &v);

 protected:
   std::vector<vm::Vec3<double>> speed_;
   size_t fiber_idx_;

   // const double k_max_speed_ = 0.1; // TODO: should be dependend on min obj
   // size
};
} // namespace geometry

#endif // FIBER_CLASS_HPP_
