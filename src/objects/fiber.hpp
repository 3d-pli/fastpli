#ifndef SRC_OBJECTS_FIBER_HPP_
#define SRC_OBJECTS_FIBER_HPP_

#include <cassert>
#include <vector>

#include "include/aabb.hpp"
#include "include/vemath.hpp"

namespace object {

class Fiber {
 public:
   explicit Fiber(const std::vector<double> &data);
   Fiber(const std::vector<vm::Vec3<double>> &points,
         const std::vector<double> &radii);
   Fiber(const std::vector<double> &points, const std::vector<double> &radii);

   // defaults
   Fiber() = default;
   Fiber(Fiber &&) = default;
   Fiber(const Fiber &) = default;
   Fiber &operator=(Fiber &&) = default;
   Fiber &operator=(const Fiber &) = default;

   ~Fiber() = default;

   // getter
   const std::vector<vm::Vec3<double>> &points() const { return points_; }
   const std::vector<double> &radii() const { return radii_; }
   std::vector<double> vector() const;
   const aabb::AABB<double, 3> &aabb() const { return aabb_; }

   size_t size() const { return points_.size(); }
   double CalcRadius(size_t idx, double t) const;

   // manipolator
   void Rotate(const vm::Mat3x3<double> &rot_mat);
   void RotateAroundPoint(const vm::Mat3x3<double> &rot_mat,
                          const vm::Vec3<double> &point);
   void Translate(const vm::Vec3<double> &translation);

   void Resize(const double f);
   void ResizePoints(const double f);
   void ResizeRadii(const double f);

 protected:
   std::vector<vm::Vec3<double>> points_;
   std::vector<double> radii_;
   aabb::AABB<double, 3> aabb_;

   void CalculateVoi();
};

using FiberBundle = std::vector<Fiber>;
using FiberBundles = std::vector<FiberBundle>;

} // namespace object

#endif // SRC_OBJECTS_FIBER_HPP_
