#ifndef FIBER_HPP_
#define FIBER_HPP_

#include <cassert>
#include <vector>

#include "include/aabb.hpp"
#include "include/vemath.hpp"

namespace data { // TODO: fiber::data

class Fiber {
 public:
   Fiber(const std::vector<vm::Vec3<float>> &points,
         const std::vector<float> &radii);
   Fiber(const std::vector<float> &points, const std::vector<float> &radii);

   // defaults
   Fiber() = default;
   Fiber(Fiber &&) = default;
   Fiber(const Fiber &) = default;
   Fiber &operator=(Fiber &&) = default;
   Fiber &operator=(const Fiber &) = default;

   ~Fiber() = default;

   // getter
   const std::vector<vm::Vec3<float>> &points() const { return points_; }
   const std::vector<float> &radii() const { return radii_; }
   const aabb::AABB<float, 3> &voi() const { return voi_; }

   size_t size() const { return points_.size(); };
   float CalcRadius(size_t idx, float t) const;

   // manipolator
   void Rotate(const std::array<float, 9> &rot_mat) {
      Rotate(vm::Mat3x3<float>(rot_mat));
   };
   void Rotate(const vm::Mat3x3<float> &rot_mat);

   void RotateAroundPoint(const std::array<float, 9> &rot_mat,
                          const std::array<float, 3> &point) {
      RotateAroundPoint(vm::Mat3x3<float>(rot_mat), vm::Vec3<float>(point));
   };

   void RotateAroundPoint(const vm::Mat3x3<float> &rot_mat,
                          const vm::Vec3<float> &point);

   void Translate(const std::array<float, 3> &translation) {
      Translate(vm::Vec3<float>(translation));
   };

   void Translate(const vm::Vec3<float> &translation);

   void Resize(const float f);
   void ResizePoints(const float f);
   void ResizeRadii(const float f);

 protected:
   std::vector<vm::Vec3<float>> points_;
   std::vector<float> radii_;
   aabb::AABB<float, 3> voi_;

   void CalculateVoi();
};

using FiberBundle = std::vector<Fiber>;
using FiberTissue = std::vector<FiberBundle>;

} // namespace data

#endif // FIBER_HPP_
