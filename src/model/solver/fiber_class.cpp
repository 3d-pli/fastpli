#include "fiber_class.hpp"

#include <cassert>
#include <vector>

#include "include/vemath.hpp"
#include "objects/fiber.hpp"

namespace object {

Fiber::Fiber(const data::Fiber &fiber_data, const size_t f_idx)
    : data::Fiber(fiber_data) {

   fiber_idx_ = f_idx;
   speed_.assign(points_.size(), vm::Vec3<float>(0));
}

size_t Fiber::ConeSize() const {
   return points_.size() < 2 ? 0 : points_.size() - 1;
}

object::Cone Fiber::Cone(size_t i) const {
   assert(i < ConeSize());
   return object::Cone(this->points_[i], points_[i + 1], radii_[i],
                       radii_[i + 1], fiber_idx_, i);
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
   for (size_t i = 0; i < points_.size(); i++) {
      auto const norm = vm::length(speed_[i]);
      if (norm > k_max_speed_)
         speed_[i] *= k_max_speed_ / norm;

      points_[i] += speed_[i];
      speed_[i] *= drag;
   }
}

bool Fiber::CheckRadius(const float obj_min_radius) {
   auto solved = true;

   if (points_.size() <= 2)
      return solved;

   if (obj_min_radius == 0)
      return solved;

   auto pos_new = points_;

   for (size_t i = 1; i < points_.size() - 1; i++) {
      auto const &p1 = points_[i - 1];
      auto const &p2 = points_[i];
      auto const &p3 = points_[i + 1];

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

         pos_new[i] = points_[i] + dv * k_max_speed_ * 0.1;

         solved = false;
      }
   }
   points_ = std::move(pos_new);

   return solved;
}

bool Fiber::CheckLength(const float obj_mean_length) {
   auto solved = true;

   if (obj_mean_length == 0)
      return solved;

   auto const min = 2.0f / 3.0f * obj_mean_length;
   auto const max = 4.0f / 3.0f * obj_mean_length;

   for (size_t i = 0; i < points_.size() - 1; i++) {
      auto distance = vm::length(points_[i + 1] - points_[i]);
      if (distance > max) {
         Split(i);
         solved = false;
      } else if (distance < min) {
         if (points_.size() == 2)
            return true;
         Combine(i);
         solved = false;
      }
   }

   return solved;
}

void Fiber::Split(size_t idx) {
   auto const pos_new = (points_[idx] + points_[idx + 1]) * 0.5;
   auto const r_new = (radii_[idx] + radii_[idx + 1]) * 0.5;
   auto const v_new = (speed_[idx] + speed_[idx + 1]) * 0.5;

   points_.insert(points_.begin() + idx + 1, pos_new);
   radii_.insert(radii_.begin() + idx + 1, r_new);
   speed_.insert(speed_.begin() + idx + 1, v_new);
}

void Fiber::Combine(size_t idx) {
   if (ConeSize() <= 1)
      return;

   if (idx == 0) {
      // don't erase first point
      points_.erase(points_.begin() + 1);
      radii_.erase(radii_.begin() + 1);
      speed_.erase(speed_.begin() + 1);
   } else if (idx == points_.size() - 2) {
      // don't erase last point
      points_.erase(points_.end() - 2);
      radii_.erase(radii_.end() - 2);
      speed_.erase(speed_.end() - 2);
   } else {
      // erase midpoint
      points_.erase(points_.begin() + idx + 1);
      radii_.erase(radii_.begin() + idx + 1);
      speed_.erase(speed_.begin() + idx + 1);

      // TODO: check
      // auto const pos_new = (points_[idx] + points_[idx + 1]) * 0.5;
      // auto const r_new = (radii_[idx] + radii_[idx + 1]) * 0.5;
      // auto const v_new = (speed_[idx] + speed_[idx + 1]) * 0.5;

      // points_.erase(points_.begin() + idx, points_.begin() + idx + 2);
      // radii_.erase(radii_.begin() + idx, radii_.begin() + idx + 2);
      // speed_.erase(speed_.begin() + idx, speed_.begin() + idx + 2);

      // points_.insert(points_.begin() + idx, pos_new);
      // radii_.insert(radii_.begin() + idx, r_new);
      // speed_.insert(speed_.begin() + idx, v_new);
   }
}

void Fiber::AddSpeed(size_t idx, const vm::Vec3<float> &v) { speed_[idx] += v; }

} // namespace object
