#include "cone_class.hpp"

#include <random>
#include <tuple>
#include <vector>

#include "include/aabb.hpp"
#include "include/vemath.hpp"

namespace object {
aabb::AABB<float, 3> Cone::aabb() const {
   return aabb::AABB<float, 3>(p0 - r, p1 + r);
}

bool Cone::CollideWith(const Cone cone) const {
   // TODO: radius varies along the length
   vm::Vec3<float> P, Q;
   std::tie(P, Q) = MinDistanceVectorCones(cone);
   float min_dist = vm::length(P - Q);
   if (min_dist < (r + cone.r))
      return true;

   return false;
}

std::tuple<vm::Vec3<float>, vm::Vec3<float>>
Cone::MinDistanceVectorCones(const Cone cone) const {
   // https://www.john.geek.nz/2009/03/code-shortest-distance-between-any-two-line-segments/
   auto const &P0 = p0;
   auto const &P1 = p1;
   auto const &Q0 = cone.p0;
   auto const &Q1 = cone.p1;

   auto const u = P1 - P0;
   auto const v = Q1 - Q0;
   auto const w = P0 - Q0;

   auto const a = vm::dot(u, u);
   auto const b = vm::dot(u, v);
   auto const c = vm::dot(v, v);
   auto const d = vm::dot(u, w);
   auto const e = vm::dot(v, w);
   auto const f = a * c - b * b;

   float sc, sN, sD = f;
   float tc, tN, tD = f;

   if (f < 1e-5f) {
      sN = 0.0;
      sD = 1.0;
      tN = e;
      tD = c;
   } else {
      sN = b * e - c * d;
      tN = a * e - b * d;
      if (sN < 0.0) {
         sN = 0.0;
         tN = e;
         tD = c;
      } else if (sN > sD) {
         sN = sD;
         tN = e + b;
         tD = c;
      }
   }

   if (tN < 0.0) {
      tN = 0.0;

      if (-d < 0.0)
         sN = 0.0;
      else if (-d > a)
         sN = sD;
      else {
         sN = -d;
         sD = a;
      }
   } else if (tN > tD) {
      tN = tD;
      if ((-d + b) < 0.0)
         sN = 0;
      else if ((-d + b) > a)
         sN = sD;
      else {
         sN = (-d + b);
         sD = a;
      }
   }

   if (std::abs(sN) < 1e-5f)
      sc = 0.0;
   else
      sc = sN / sD;
   if (std::abs(tN) < 1e-5f)
      tc = 0.0;
   else
      tc = tN / tD;

   // vm::Vec3<float> dP;
   // dP = w + (u * sc) - (v * tc);
   auto P = P0 + u * sc;
   auto Q = Q0 + v * tc;

   return std::make_tuple(P, Q);
}

std::tuple<vm::Vec3<float>, vm::Vec3<float>, vm::Vec3<float>, vm::Vec3<float>>
Cone::PushConesApart(const Cone cone) const {

   // TODO: clean up

   vm::Vec3<float> P, Q;
   std::tie(P, Q) = MinDistanceVectorCones(cone);

   auto delta = P - Q;
   float norm = vm::length(P - Q);

   if (norm < 1e-8f) {
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_real_distribution<float> dis(-1.0, 1.0);
      delta = vm::Vec3<float>(dis(gen), dis(gen), dis(gen));

      return std::make_tuple(
          vm::Vec3<float>(dis(gen), dis(gen), dis(gen)) * 0.1,
          vm::Vec3<float>(dis(gen), dis(gen), dis(gen)) * 0.1,
          vm::Vec3<float>(dis(gen), dis(gen), dis(gen)) * 0.1,
          vm::Vec3<float>(dis(gen), dis(gen), dis(gen)) * 0.1);
   }

   auto delta_speed =
       0.05 * std::min(r0, std::min(r1, std::min(cone.r0, cone.r1)));

   // auto f = (P - Q) / vm::length(P - Q) * delta_speed;

   auto f0 =
       delta / norm * delta_speed * vm::length(P - p1) / vm::length(p1 - p0);
   auto f1 =
       delta / norm * delta_speed * vm::length(P - p0) / vm::length(p1 - p0);

   auto f2 = -delta / norm * delta_speed * vm::length(Q - cone.p1) /
             vm::length(cone.p1 - cone.p0);
   auto f3 = -delta / norm * delta_speed * vm::length(Q - cone.p0) /
             vm::length(cone.p1 - cone.p0);

   return std::make_tuple(f0, f1, f2, f3);
   // return std::make_tuple(f, f, -f, -f);
}

} // namespace object
