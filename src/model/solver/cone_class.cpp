#include "cone_class.hpp"

#include <random>
#include <tuple>
#include <vector>

#include "aabb.hpp"
#include "vemath.hpp"

namespace object {
aabb::AABB<float, 3> Cone::aabb() const {
   return aabb::AABB<float, 3>(p0 - r, p1 + r);
}

bool Cone::CollideWith(const Cone cone) const {
   // TODO: radius varies along the length
   float min_dist = MinDistanceCones(cone);
   if (min_dist < (r + cone.r))
      return true;

   return false;
}

float Cone::MinDistanceCones(const Cone cone) const {
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

   if (f < 1e-9f) {
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

   if (std::abs(sN) < 1e-9f)
      sc = 0.0;
   else
      sc = sN / sD;
   if (std::abs(tN) < 1e-9f)
      tc = 0.0;
   else
      tc = tN / tD;

   auto dP = w + (u * sc) - (v * tc);
   return vm::length(dP);
}

vm::Vec3<float> Cone::MinDistanceVectorCones(const Cone cone) const {
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

   if (f < 1e-9f) {
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

   if (std::abs(sN) < 1e-9f)
      sc = 0.0;
   else
      sc = sN / sD;
   if (std::abs(tN) < 1e-9f)
      tc = 0.0;
   else
      tc = tN / tD;

   vm::Vec3<float> dP;
   dP = w + (u * sc) - (v * tc);

   return dP;
}

vm::Vec3<float> Cone::PushConesApart(const Cone cone) const {

   auto delta = MinDistanceVectorCones(cone);
   float norm = vm::length(delta);

   if (norm < 1e-8f) {
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_real_distribution<float> dis(-1.0, 1.0);
      delta = vm::Vec3<float>(dis(gen), dis(gen), dis(gen));
   }

   vm::normalize(delta);
   delta *= 0.5;

   return delta * 0.1;
}

} // namespace object
