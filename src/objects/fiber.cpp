#include "fiber.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

#include "include/aabb.hpp"
#include "include/vemath.hpp"

namespace object {
Fiber::Fiber(const std::vector<double> &data) {

   if (data.size() % 4 != 0)
      throw std::invalid_argument("fiber data has to be of shape nx4");

   points_.resize(data.size() / 4);
   points_.shrink_to_fit();
   radii_.resize(data.size() / 4);
   radii_.shrink_to_fit();
   aabb_ = aabb::AABB<double, 3>{};

   for (size_t i = 0; i < data.size() / 4; i++) {
      points_[i] =
          vm::Vec3<double>(data[4 * i + 0], data[4 * i + 1], data[4 * i + 2]);
      radii_[i] = data[4 * i + 3];
   }

   // calc voi
   if (points_.empty()) {
      return;
   } else if (points_.size() == 1) {
      aabb_ = aabb::AABB<double, 3>(points_[0], points_[0]);
      return;
   }
   CalculateVoi();
}

Fiber::Fiber(const std::vector<double> &points,
             const std::vector<double> &radii) {

   if (points.size() != radii.size() * 3)
      throw std::invalid_argument("points and radii aren't same size");

   points_.resize(radii.size());
   points_.shrink_to_fit();
   aabb_ = aabb::AABB<double, 3>{};
   radii_ = radii;

   for (size_t i = 0; i < radii_.size(); i++)
      points_[i] = vm::Vec3<double>(points[3 * i + 0], points[3 * i + 1],
                                    points[3 * i + 2]);

   // calc voi
   if (points_.empty()) {
      return;
   } else if (points_.size() == 1) {
      aabb_ = aabb::AABB<double, 3>(points_[0], points_[0]);
      return;
   }
   CalculateVoi();
}

Fiber::Fiber(const std::vector<vm::Vec3<double>> &points,
             const std::vector<double> &radii) {

   if (points.size() != radii.size())
      throw std::invalid_argument("points and radii aren't same size");

   points_ = points;
   radii_ = radii;

   // calc voi
   if (points_.empty()) {
      return;
   } else if (points_.size() == 1) {
      aabb_ = aabb::AABB<double, 3>(points_[0], points_[0]);
      return;
   }
   CalculateVoi();
}

void Fiber::CalculateVoi() {

   if (points_.empty()) {
      aabb_ = aabb::AABB<double, 3>{};
      return;
   }

   aabb_ =
       aabb::AABB<double, 3>(points_[0] - radii_[0], points_[0] + radii_[0]);
   for (size_t i = 0; i < points_.size() - 1; i++) {
      auto r_max = std::max(radii_[i], radii_[i + 1]);
      auto tmp_voi = aabb::AABB<double, 3>(points_[i], points_[i + 1]);
      tmp_voi.min -= r_max;
      tmp_voi.max += r_max;
      aabb_.Unite(tmp_voi);
   }
}

std::vector<double> Fiber::vector() const {
   std::vector<double> fiber;
   fiber.reserve(points_.size() * 4);

   for (size_t i = 0; i < points_.size(); i++) {
      fiber.push_back(points_[i].x());
      fiber.push_back(points_[i].y());
      fiber.push_back(points_[i].z());
      fiber.push_back(radii_[i]);
   }

   return fiber;
}

double Fiber::CalcRadius(size_t idx, double t) const {
   assert(t >= 0 && t <= 1);
   assert(idx + 1 < radii_.size());

   return radii_[idx] * (1 - t) + radii_[idx + 1] * t;
}

void Fiber::Rotate(const vm::Mat3x3<double> &rot_mat) {
   for (auto &p : points_)
      p = vm::dot(rot_mat, p);
   CalculateVoi();
}

void Fiber::Translate(const vm::Vec3<double> &translation) {
   auto off = vm::Vec3<double>(translation);
   for (auto &p : points_)
      p += off;
   CalculateVoi();
}

void Fiber::RotateAroundPoint(const vm::Mat3x3<double> &rot_mat,
                              const vm::Vec3<double> &point) {
   auto rot = vm::Mat3x3<double>(rot_mat);
   auto off = vm::Vec3<double>(point);
   for (auto &p : points_)
      p = vm::dot(rot, p - off) + off;
   CalculateVoi();
}

void Fiber::Resize(const double f) {
   for (auto &p : points_)
      p *= f;
   for (auto &r : radii_)
      r *= f;
   CalculateVoi();
}
void Fiber::ResizePoints(const double f) {
   for (auto &p : points_)
      p *= f;
   CalculateVoi();
}

void Fiber::ResizeRadii(const double f) {
   for (auto &r : radii_)
      r *= f;
   CalculateVoi();
}

} // namespace object
