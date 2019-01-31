#include "fiber_bundle.hpp"

#include <cassert>
#include <iostream>
#include <vector>

#include "include/aabb.hpp"
#include "include/vemath.hpp"
#include "objects/fiber.hpp"

namespace fiber {

namespace layer {

Property::Property(float s, float n, float m, char o) : scale(s), dn(n), mu(m) {
   if (o == 'b')
      orientation = Orientation::background;
   else if (o == 'p')
      orientation = Orientation::parallel;
   else if (o == 'r')
      orientation = Orientation::radial;
   else
      throw std::invalid_argument("Orientation must be \"b\", \"p\", or \"r\"");
}

size_t Properties::size() const {
   assert(scale_.size() == scale_sqr_.size());
   assert(scale_.size() == dn_.size());
   assert(scale_.size() == mu_.size());
   assert(scale_.size() == orientation_.size());
   return scale_.size();
}

void Properties::clear() {
   scale_.clear();
   scale_sqr_.clear();
   dn_.clear();
   mu_.clear();
   orientation_.clear();
}

void Properties::reserve(size_t i) {
   scale_.reserve(i);
   scale_sqr_.reserve(i);
   dn_.reserve(i);
   mu_.reserve(i);
   orientation_.reserve(i);
};

void Properties::resize(size_t i) {
   scale_.resize(i);
   scale_sqr_.resize(i);
   dn_.resize(i);
   mu_.resize(i);
   orientation_.resize(i);
};

void Properties::push_back(Property p) {
   scale_.push_back(p.scale);
   scale_sqr_.push_back(p.scale * p.scale);
   dn_.push_back(p.dn);
   mu_.push_back(p.mu);
   orientation_.push_back(p.orientation);
}
} // namespace layer

/**
 * FiberBundle
 */
Bundle::Bundle(std::vector<Fiber> fibers,
               std::vector<layer::Property> properties) {
   fibers_.swap(fibers);

   // sort for layer size
   std::sort(properties.begin(), properties.end(),
             [&](layer::Property p1, layer::Property p2) {
                return p1.scale < p2.scale;
             });

   layers_.reserve(properties.size());
   for (auto const &p : properties)
      layers_.push_back(p);
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

void Bundle::Rotate(const vm::Mat3x3<float> &rot_mat) {
   for (auto &fiber : fibers_)
      fiber.Rotate(rot_mat);
}

void Bundle::RotateAroundPoint(const vm::Mat3x3<float> &rot_mat,
                               const vm::Vec3<float> &point) {
   for (auto &fiber : fibers_)
      fiber.RotateAroundPoint(rot_mat, point);
}

void Bundle::Translate(const vm::Vec3<float> &translation) {
   for (auto &fiber : fibers_)
      fiber.Translate(translation);
}

} // namespace fiber
