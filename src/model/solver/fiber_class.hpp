#ifndef FIBER_CLASS_HPP_
#define FIBER_CLASS_HPP_

#include <vector>

#include "fiber_segment.hpp"
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

   // fiber segments
   size_t FiberSegmentSize() const;
   geometry::FiberSegment FiberSegment(size_t i) const;
   std::vector<geometry::FiberSegment> FiberSegments() const;

   // manipulators
   const double &max_speed() const { return max_speed_; };
   void set_max_speed(double max_speed) { max_speed_ = max_speed; };

   void Move();
   void Drag(const double drag = 1);
   bool ApplyCurvatureConstrain(const double obj_min_radius);
   bool ApplyFiberSegmentLengthConstrain(const double obj_mean_length);

   void Split(size_t idx);
   void Combine(size_t idx);

   void AddSpeed(size_t idx, const vm::Vec3<double> &v);

 protected:
   size_t fiber_idx_;
   double max_speed_ = 0;
   std::vector<vm::Vec3<double>> speed_;
};
} // namespace geometry

#endif // FIBER_CLASS_HPP_
