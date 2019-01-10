#ifndef FIBER_CLASS_HPP_
#define FIBER_CLASS_HPP_

#include <vector>

#include "cone_class.hpp"
#include "include/vemath.hpp"
#include "objects/fiber.hpp"

namespace geometry {
class Fiber : public object::fiber::Fiber {
 public:
   // defaults
   Fiber(Fiber &&) = default;
   Fiber(const Fiber &) = default;
   Fiber &operator=(Fiber &&) = default;
   Fiber &operator=(const Fiber &) = default;
   ~Fiber() = default;

   // constructors
   Fiber(const object::fiber::Fiber &fiber, const size_t fiber_idx);

   // getter
   const size_t &fiber_idx() const { return fiber_idx_; }
   const std::vector<vm::Vec3<float>> &speed() const { return speed_; }

   // cones
   size_t ConeSize() const;
   object::Cone Cone(size_t i) const;
   std::vector<object::Cone> Cones() const;

   // manipulators
   void Move(const float drag = 1);
   bool CheckRadius(const float obj_min_radius);
   bool CheckLength(const float obj_mean_length);

   void Split(size_t idx);
   void Combine(size_t idx);

   void AddSpeed(size_t idx, const vm::Vec3<float> &v);

 protected:
   std::vector<vm::Vec3<float>> speed_;
   size_t fiber_idx_;
   const float k_max_speed_ = 0.1; // TODO: should be dependend on min obj size
};
} // namespace geometry

#endif // FIBER_CLASS_HPP_
