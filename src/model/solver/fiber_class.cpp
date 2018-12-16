#include "fiber_class.hpp"

#include <cassert>
#include <vector>

#include "vemath.hpp"

namespace geometry {

Fiber::Fiber(const std::vector<vm::Vec3<float>> &p, const std::vector<float> &r,
             const size_t f_idx)
    : FiberRawData<float>(p, r) {

   fiber_idx = f_idx;
   speed.assign(points.size(), vm::Vec3<float>(0));
}

Fiber::Fiber(const std::vector<float> &p, const std::vector<float> &r,
             const size_t f_idx)
    : FiberRawData<float>(p, r) {

   fiber_idx = f_idx;
   speed.assign(points.size(), vm::Vec3<float>(0));
}

Fiber::Fiber(const FiberRawData<float> &fiber_data, const size_t f_idx)
    : FiberRawData<float>(fiber_data) {

   fiber_idx = f_idx;
   speed.assign(points.size(), vm::Vec3<float>(0));
}

size_t Fiber::ConeSize() const {
   return points.size() < 2 ? 0 : points.size() - 1;
}

object::Cone Fiber::Cone(size_t i) const {
   assert(i < ConeSize());
   return object::Cone(points[i], points[i + 1], radii[i], radii[i + 1],
                       fiber_idx, i);
}

std::vector<object::Cone> Fiber::Cones() const {
   std::vector<object::Cone> data;
   data.reserve(ConeSize());
   for (size_t i = 0; i < ConeSize(); i++) {
      data.push_back(Cone(i));
   }
   return data;
}

void Fiber::Move(const float drag) {
   for (size_t i = 0; i < points.size(); i++) {
      auto const norm = vm::length(speed[i]);
      if (norm > k_max_speed_)
         speed[i] *= k_max_speed_ / norm;

      points[i] += speed[i];
      speed[i] *= drag;
   }
}

bool Fiber::CheckRadius(const float obj_min_radius) {
   auto solved = true;

   if (points.size() <= 2)
      return solved;

   if (obj_min_radius == 0)
      return solved;

   auto pos_new = points;

   for (size_t i = 1; i < points.size() - 1; i++) {
      auto const &p1 = points[i - 1];
      auto const &p2 = points[i];
      auto const &p3 = points[i + 1];

      auto const a = vm::length(p1 - p2);
      auto const b = vm::length(p2 - p3);
      auto const c = vm::length(p3 - p1);

      auto const s = (a + b + c) / 2.0f;
      auto const f = std::sqrt(s * (s - a) * (s - b) * (s - c));
      auto const r = a * b * c / (4 * f);

      // TODO: min_radius sollte gr;-er sein als gleichseitiges dreieck,
      // R=sqrt(3)/3.0*mean_...
      // TODO: therefore if angle < 60 dann 'ndern
      auto const gamma = std::acos((c * c - a * a - b * b) / (-2.0f * a * b));

      if (r < obj_min_radius || gamma < 60.0f / 180.0f * M_PI) {
         // std::cout << "i:" << i << " " << r << "<" << obj_min_radius << ", "
         // << k_max_speed_ << std::endl;
         auto const v1 = p1 - p2;
         auto const v2 = p3 - p2;
         auto dv = v1 + v2;
         vm::normalize(dv);

         pos_new[i] = points[i] + dv * k_max_speed_ * 0.1;

         solved = false;
      }
   }
   points = std::move(pos_new);

   return solved;
}

bool Fiber::CheckLength(const float obj_mean_length) {
   auto solved = true;

   if (points.size() <= 2)
      return solved;

   if (obj_mean_length == 0)
      return solved;

   auto const min = 2.0f / 3.0f * obj_mean_length;
   auto const max = 4.0f / 3.0f * obj_mean_length;

   for (size_t i = 0; i < points.size() - 1; i++) {
      auto distance = vm::length(points[i + 1] - points[i]);
      if (distance > max) {
         Split(i);
         solved = false;
      } else if (distance < min) {
         if (points.size() == 2)
            return true;
         Combine(i);
         solved = false;
      }
   }

   return solved;
}

void Fiber::Split(size_t idx) {
   auto const pos_new = (points[idx] + points[idx + 1]) * 0.5;
   auto const r_new = (radii[idx] + radii[idx + 1]) * 0.5;
   auto const v_new = (speed[idx] + speed[idx + 1]) * 0.5;

   points.insert(points.begin() + idx + 1, pos_new);
   radii.insert(radii.begin() + idx + 1, r_new);
   speed.insert(speed.begin() + idx + 1, v_new);
}

void Fiber::Combine(size_t idx) {
   if (ConeSize() <= 1)
      return;

   if (idx == 0) {
      // don't erase first point
      points.erase(points.begin() + 1);
      radii.erase(radii.begin() + 1);
      speed.erase(speed.begin() + 1);
   } else if (idx == points.size() - 2) {
      // don't erase last point
      points.erase(points.end() - 2);
      radii.erase(radii.end() - 2);
      speed.erase(speed.end() - 2);
   } else {
      // erase midpoint
      points.erase(points.begin() + idx + 1);
      radii.erase(radii.begin() + idx + 1);
      speed.erase(speed.begin() + idx + 1);

      // TODO: check
      // auto const pos_new = (points[idx] + points[idx + 1]) * 0.5;
      // auto const r_new = (radii[idx] + radii[idx + 1]) * 0.5;
      // auto const v_new = (speed[idx] + speed[idx + 1]) * 0.5;

      // points.erase(points.begin() + idx, points.begin() + idx + 2);
      // radii.erase(radii.begin() + idx, radii.begin() + idx + 2);
      // speed.erase(speed.begin() + idx, speed.begin() + idx + 2);

      // points.insert(points.begin() + idx, pos_new);
      // radii.insert(radii.begin() + idx, r_new);
      // speed.insert(speed.begin() + idx, v_new);
   }
}

} // namespace object
