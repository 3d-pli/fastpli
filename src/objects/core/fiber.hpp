#ifndef FIBER_HPP_
#define FIBER_HPP_

#include <cassert>
#include <vector>

#include "aabb.hpp"
#include "vemath.hpp"

namespace object {

class FiberData {
 public:
   FiberData(const std::vector<vm::Vec3<float>> &points,
             const std::vector<float> &radii);
   FiberData(const std::vector<float> &points, const std::vector<float> &radii);

   // defaults
   FiberData() = default;
   FiberData(FiberData &&) = default;
   FiberData(const FiberData &) = default;
   FiberData &operator=(FiberData &&) = default;
   FiberData &operator=(const FiberData &) = default;

   ~FiberData() = default;

   // getter
   // TODO: return of arrray should be only in python module
   const std::vector<vm::Vec3<float>> &points() const { return points_; }
   const std::vector<float> &radii() const { return radii_; }
   const aabb::AABB<float, 3> &voi() const { return voi_; }

   size_t size() const { return points_.size(); };
   float CalcRadius(size_t idx, float t) const;

   // manipolator
   void RotateFiber(const std::array<float, 9> &rot_mat) {
      RotateFiber(vm::Mat3x3<float>(rot_mat));
   };
   void RotateFiber(const vm::Mat3x3<float> &rot_mat);

   void RotateFiberAroundPoint(const std::array<float, 9> &rot_mat,
                               const std::array<float, 3> &point) {
      RotateFiberAroundPoint(vm::Mat3x3<float>(rot_mat),
                             vm::Vec3<float>(point));
   };

   void RotateFiberAroundPoint(const vm::Mat3x3<float> &rot_mat,
                               const vm::Vec3<float> point);

   void TranslateFiber(const std::array<float, 3> &translation) {
      TranslateFiber(vm::Vec3<float>(translation));
   };

   void TranslateFiber(const vm::Vec3<float> &translation);

   void ResizeFiberPos(const float f);
   void ResizeFiberRadii(const float f);
   void ResizeFiber(const float f);

 private:
   std::vector<vm::Vec3<float>> points_;
   std::vector<float> radii_;
   aabb::AABB<float, 3> voi_;

   void CalculateVoi();
};

} // namespace object

#endif // FIBER_HPP_
