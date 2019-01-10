#ifndef SIMULATION_FIBER_MODEL_HPP_
#define SIMULATION_FIBER_MODEL_HPP_

#include <utility>
#include <vector>

#include "include/aabb.hpp"
#include "include/vemath.hpp"
#include "objects/fiber.hpp"

namespace fiber {

using Data = data::Fiber;

namespace layer {
enum class Orientation { background, parallel, radial };

struct Property {
   Property(float s, float n, float m, char o) : scale(s), dn(n), mu(m) {
      if (o == 'b')
         orientation = Orientation::background;
      else if (o == 'p')
         orientation = Orientation::parallel;
      else if (o == 'r')
         orientation = Orientation::radial;
      else
         throw std::invalid_argument(
             "Orientation must be \"b\", \"p\", or \"r\"");
   }
   Property(float s, float n, float m, Orientation o)
       : scale(s), dn(n), mu(m), orientation(o) {}

   float scale{};
   float dn{};
   float mu{};
   Orientation orientation{Orientation::background};
};
} // namespace layer

class Bundle {
 public:
   Bundle(std::vector<Data> fibers, std::vector<layer::Property> properties);
   ~Bundle() = default;

   // getter
   const Data &fiber(size_t i) const { return fibers_[i]; }
   const std::vector<Data> &fibers() const { return fibers_; }
   size_t size() const { return fibers_.size(); }
   const std::vector<float> &layer_scale() const { return layer_scale_; }
   const std::vector<float> &layer_scale_sqr() const {
      return layer_scale_sqr_;
   }

   size_t layer_size() const { return layer_dn_.size(); }
   float layer_dn(size_t i) const { return layer_dn_[i]; }
   float layer_mu(size_t i) const { return layer_mu_[i]; }
   layer::Orientation layer_orientation(size_t i) const {
      return layer_orientation_[i];
   }
   const std::vector<float> &layer_dn() const { return layer_dn_; }
   const std::vector<float> &layer_mu() const { return layer_mu_; }
   const std::vector<layer::Orientation> &layer_orientation() const {
      return layer_orientation_;
   }
   const aabb::AABB<float, 3> &voi() const { return voi_; }

   // manipulator
   void Resize(const float f);
   void ResizePoints(const float f);
   void ResizeRadii(const float f);
   void Rotate(const std::array<float, 9> &rot_mat);
   void RotateAroundPoint(const std::array<float, 9> &rot_mat,
                          std::array<float, 3> point);
   void Translate(const std::array<float, 3> &translation);

 private:
   std::vector<Data> fibers_;
   std::vector<float> layer_scale_;
   std::vector<float> layer_scale_sqr_;
   std::vector<float> layer_dn_;
   std::vector<float> layer_mu_;
   std::vector<layer::Orientation> layer_orientation_;

   aabb::AABB<float, 3> voi_{};
};
} // namespace fiber

#endif // FIBER_MODEL_HPP_
