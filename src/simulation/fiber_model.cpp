#include "fiber_model.hpp"

#include <cassert>
#include <iostream>
#include <vector>

#include "include/aabb.hpp"
#include "include/vemath.hpp"
#include "objects/fiber.hpp"

namespace fiber {
Bundle::Bundle(std::vector<data::Fiber> fibers,
               std::vector<layer::Property> properties) {
   fibers_.swap(fibers);

   std::vector<float> layer_dn, layer_mu, layer_scale;
   std::vector<layer::Orientation> layer_orientation;

   layer_dn.reserve(properties.size());
   layer_mu.reserve(properties.size());
   layer_scale.reserve(properties.size());
   layer_orientation.reserve(properties.size());

   for (const auto &p : properties) {
      layer_dn.push_back(p.dn);
      layer_mu.push_back(p.mu);
      layer_scale.push_back(p.scale);
      layer_orientation.push_back(p.orientation);
   }

   // sort layer
   std::vector<size_t> idx(layer_scale.size());
   std::iota(idx.begin(), idx.end(), 0);
   std::sort(idx.begin(), idx.end(), [&](size_t i1, size_t i2) {
      return layer_scale[i1] < layer_scale[i2];
   });

   layer_dn_.reserve(properties.size());
   layer_mu_.reserve(properties.size());
   layer_scale_.reserve(properties.size());
   layer_scale_sqr_.reserve(properties.size());
   layer_orientation_.reserve(properties.size());

   for (auto i : idx) {
      layer_dn_.push_back(layer_dn[i]);
      layer_mu_.push_back(layer_mu[i]);
      layer_scale_.push_back(layer_scale[i]);
      layer_scale_sqr_.push_back(layer_scale[i] * layer_scale[i]);
      layer_orientation_.push_back(layer_orientation[i]);
   }
}

void Bundle::Resize(const float f) {
   for (auto &fiber : fibers_)
      fiber.Resize(f);
}
void Bundle::ResizePoints(const float f) {
   for (auto &fiber : fibers_)
      fiber.ResizePoints(f);
}

void Bundle::ResizeRadii(const float f) {
   for (auto &fiber : fibers_)
      fiber.ResizeRadii(f);
}

void Bundle::Rotate(const std::array<float, 9> &rot_mat) {
   for (auto &fiber : fibers_)
      fiber.Rotate(rot_mat);
}

void Bundle::RotateAroundPoint(const std::array<float, 9> &rot_mat,
                               const std::array<float, 3> point) {
   for (auto &fiber : fibers_)
      fiber.RotateAroundPoint(rot_mat, point);
}

void Bundle::Translate(const std::array<float, 3> &translation) {
   for (auto &fiber : fibers_)
      fiber.Translate(translation);
}

} // namespace fiber
