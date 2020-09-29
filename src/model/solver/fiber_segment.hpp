#ifndef FIBER_SEGMENT_HPP_
#define FIBER_SEGMENT_HPP_

#include <tuple>
#include <vector>

#include "include/aabb.hpp"
#include "include/vemath.hpp"

namespace geometry {
struct FiberSegment {
   // geometrically forms a conical obj

   // defaults
   FiberSegment() = default;
   FiberSegment(FiberSegment &&) = default;
   FiberSegment(const FiberSegment &) = default;
   FiberSegment &operator=(FiberSegment &&) = default;
   FiberSegment &operator=(const FiberSegment &) = default;
   ~FiberSegment() = default;

   // constructors
   FiberSegment(const vm::Vec3<double> p0, const vm::Vec3<double> p1,
                const double r0, const double r1, const long long fiber_idx,
                const long long fiber_elm)
       : p0(p0), p1(p1), r0(r0), r1(r1), r(r0 > r1 ? r0 : r1),
         fiber_idx(fiber_idx), fiber_elm(fiber_elm){};

   aabb::AABB<double, 3> aabb() const;

   // colliding solving functions
   bool CollideWith(const FiberSegment obj) const;
   std::tuple<vm::Vec3<double>, vm::Vec3<double>>
   MinDistanceVector(const FiberSegment obj) const;
   std::tuple<vm::Vec3<double>, vm::Vec3<double>, vm::Vec3<double>,
              vm::Vec3<double>, double>
   PushApart(const FiberSegment obj) const;

   vm::Vec3<double> p0, p1;
   double r0, r1, r; // TODO: check if r0, r1 are needed
   size_t fiber_idx, fiber_elm;
};
} // namespace geometry

#endif // FIBER_SEGMENT_HPP_
