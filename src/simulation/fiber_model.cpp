#include "fiber_model.hpp"

#include <cassert>
#include <iostream>
#include <vector>

#include "include/aabb.hpp"
#include "include/vemath.hpp"
#include "objects/fiber.hpp"

FiberData::FiberData(const std::vector<float> &points,
                     const std::vector<float> &radii) {

   if (points.size() != radii.size() * 3)
      throw std::invalid_argument("points and radii aren't same size");

   pos_.clear();
   voi_ = aabb::AABB<float, 3>{};

   radii_ = radii;
   pos_.reserve(radii_.size());
   for (size_t i = 0; i < radii_.size(); i += 3)
      pos_.push_back(
          vm::Vec3<double>(points[i + 0], points[i + 1], points[i + 2]));

   assert(pos_.size() == radii_.size());

   // calc voi
   if (pos_.empty())
      return;
   else if (pos_.size() == 1) {
      voi_ = aabb::AABB<double, 3>(pos_[0], pos_[0]);
      return;
   }
   CalculateVoi();
}

FiberData::FiberData(const std::vector<vm::Vec3<double>> &points,
                     const std::vector<double> &radii) {

   if (points.size() != radii.size())
      throw std::invalid_argument("points and radii aren't same size");

   pos_ = points;
   radii_ = radii;

   // calc voi
   if (pos_.empty())
      return;
   else if (pos_.size() == 1) {
      voi_ = aabb::AABB<double, 3>(pos_[0], pos_[0]);
      return;
   }
   CalculateVoi();
}

void FiberData::CalculateVoi() {

   if (pos_.empty()) {
      voi_ = aabb::AABB<double, 3>{};
      return;
   }

   voi_ = aabb::AABB<double, 3>(pos_[0], pos_[0]);
   for (size_t i = 0; i < pos_.size() - 1; i++) {
      auto r_max = std::max(radii_[i], radii_[i + 1]);
      auto tmp_voi = aabb::AABB<double, 3>(pos_[i], pos_[i + 1]);
      tmp_voi.min -= r_max;
      tmp_voi.max += r_max;
      voi_.Unite(tmp_voi);
   }
}

double FiberData::CalcRadius(size_t idx, double t) const {
   assert(t >= 0 && t <= 1);
   assert(idx + 1 < radii_.size());

   return radii_[idx] * (1 - t) + radii_[idx + 1] * t;
}

void FiberData::RotateFiber(const vm::Mat3x3<double> &rot_mat) {
   for (auto &p : pos_)
      p = vm::dot(rot_mat, p);
   CalculateVoi();
}

void FiberData::TranslateFiber(const vm::Vec3<double> &translation) {
   auto off = vm::Vec3<double>(translation);
   for (auto &p : pos_)
      p += off;
   CalculateVoi();
}

void FiberData::RotateFiberAroundPoint(const vm::Mat3x3<double> &rot_mat,
                                       const vm::Vec3<double> point) {
   auto rot = vm::Mat3x3<double>(rot_mat);
   auto off = vm::Vec3<double>(point);
   for (auto &p : pos_)
      p = vm::dot(rot, p - off) + off;
   CalculateVoi();
}

void FiberData::ResizeFiberPos(const double f) {
   for (auto &p : pos_)
      p *= f;
   CalculateVoi();
}

void FiberData::ResizeFiberRadii(const double f) {
   for (auto &r : radii_)
      r *= f;
   CalculateVoi();
}

void FiberData::ResizeFiber(const double f) {
   for (auto &p : pos_)
      p *= f;
   for (auto &r : radii_)
      r *= f;
   CalculateVoi();
}

/**
 * FiberBundle::
 */
void FiberBundle::push_fiber(FiberData &fiber) {
   // unify fiber voi with fiber_bundle voi
   voi_.Unite(fiber.voi());
   fibers_.push_back(fiber);
}

void FiberBundle::SetFiberBundleProperties(
    const std::vector<LayerProperty> &prop) {

   std::vector<double> layer_dn, layer_mu, layer_scale;
   std::vector<Orientation> layer_orientation;

   for (const auto &p : prop) {
      layer_dn.push_back(p.dn);
      layer_mu.push_back(p.mu);
      layer_scale.push_back(p.scale);
      layer_orientation.push_back(p.orientation);
   }

   // sort for scale
   std::vector<size_t> idx(layer_scale.size());
   std::iota(idx.begin(), idx.end(), 0);
   std::sort(idx.begin(), idx.end(), [&](size_t i1, size_t i2) {
      return layer_scale[i1] < layer_scale[i2];
   });

   for (auto i : idx) {
      layer_dn_.push_back(layer_dn[i]);
      layer_mu_.push_back(layer_mu[i]);
      layer_scale_.push_back(layer_scale[i]);
      layer_scale_sqr_.push_back(layer_scale[i] * layer_scale[i]);
      layer_orientation_.push_back(layer_orientation[i]);
   }
}

void FiberBundle::RotateFiberBundle(const std::array<double, 9> &rot_mat) {
   for (auto &fiber : fibers_)
      fiber.RotateFiber(rot_mat);
}

void FiberBundle::RotateFiberBundleAroundPoint(
    const std::array<double, 9> &rot_mat, const std::array<double, 3> point) {
   for (auto &fiber : fibers_)
      fiber.RotateFiberAroundPoint(rot_mat, point);
}

void FiberBundle::TranslateFiberBundle(
    const std::array<double, 3> &translation) {
   for (auto &fiber : fibers_)
      fiber.TranslateFiber(translation);
}

void FiberBundle::ResizeFiberBundlePos(const double f) {
   for (auto &fiber : fibers_)
      fiber.ResizeFiberPos(f);
}

void FiberBundle::ResizeFiberBundleRadii(const double f) {
   for (auto &fiber : fibers_)
      fiber.ResizeFiberRadii(f);
}

void FiberBundle::ResizeFiberBundle(const double f) {
   for (auto &fiber : fibers_)
      fiber.ResizeFiber(f);
}
