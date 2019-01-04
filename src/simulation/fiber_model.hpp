#ifndef FIBER_MODEL_HPP_
#define FIBER_MODEL_HPP_

#include <utility>
#include <vector>

#include "include/aabb.hpp"
#include "include/vemath.hpp"
#include "objects/fiber.hpp"

namespace fiber {
namespace layer {
enum class Orientation { background, parallel, radial };

struct Property {
   Property(double s, double n, double m, char o) : scale(s), dn(n), mu(m) {
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
   Property(double s, double n, double m, Orientation o)
       : scale(s), dn(n), mu(m), orientation(o) {}

   double scale{};
   double dn{};
   double mu{};
   Orientation orientation{Orientation::background};
};
} // namespace layer

class Bundle {
 public:
   Bundle(std::vector<data::Fiber> fibers,
          std::vector<layer::Property> properties);
   ~Bundle() = default;

   // getter
   const data::Fiber &fiber(size_t i) const { return fibers_[i]; }
   const std::vector<data::Fiber> &fibers() const { return fibers_; }
   size_t size() const { return fibers_.size(); }
   const std::vector<double> &layer_scale() const { return layer_scale_; }
   const std::vector<double> &layer_scale_sqr() const {
      return layer_scale_sqr_;
   }

   size_t layer_size() const { return layer_dn_.size(); }
   double layer_dn(size_t i) const { return layer_dn_[i]; }
   double layer_mu(size_t i) const { return layer_mu_[i]; }
   double layer_orientation(size_t i) const { return layer_orientation_[i]; }
   const std::vector<double> &layer_dn() const { return layer_dn_; }
   const std::vector<double> &layer_mu() const { return layer_mu_; }
   const std::vector<layer::Orientation> &layer_orientation() const {
      return layer_orientation_;
   }
   const aabb::AABB<double, 3> &voi() const { return voi_; }

 private:
   std::vector<data::Fiber> fibers_;
   std::vector<double> layer_scale_;
   std::vector<double> layer_scale_sqr_;
   std::vector<double> layer_dn_;
   std::vector<double> layer_mu_;
   std::vector<layer::Orientation> layer_orientation_;

   aabb::AABB<double, 3> voi_{};
};
} // namespace fiber

#endif // FIBER_MODEL_HPP_