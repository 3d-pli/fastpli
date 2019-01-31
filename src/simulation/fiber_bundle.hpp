#ifndef SIMULATION_FIBER_BUNDLE_HPP_
#define SIMULATION_FIBER_BUNDLE_HPP_

#include <utility>
#include <vector>

#include "include/aabb.hpp"
#include "include/vemath.hpp"
#include "objects/fiber.hpp"

namespace fiber {

namespace layer {
enum class Orientation { background, parallel, radial };

struct Property {
   Property(float s, float n, float m, char o);
   Property(float s, float n, float m, Orientation o)
       : scale(s), dn(n), mu(m), orientation(o) {}

   float scale{};
   float dn{};
   float mu{};
   Orientation orientation{};
};

class Properties {
 public:
   const std::vector<float> &scale() const { return scale_; };
   const std::vector<float> &scale_sqr() const { return scale_sqr_; };
   const std::vector<float> &dn() const { return dn_; };
   const std::vector<float> &mu() const { return mu_; };
   const std::vector<Orientation> &orientation() const { return orientation_; };

   size_t size() const;
   void clear();
   void reserve(size_t i);
   void resize(size_t i);
   void push_back(Property p);
   void push_back(float s, float n, float m, char o) {
      push_back(Property(s, n, m, o));
   };
   void push_back(float s, float n, float m, Orientation o) {
      push_back(Property(s, n, m, o));
   };

 private:
   std::vector<float> scale_;
   std::vector<float> scale_sqr_;
   std::vector<float> dn_;
   std::vector<float> mu_;
   std::vector<Orientation> orientation_;
};

} // namespace layer

using Fiber = object::fiber::Fiber;

class Bundle {
 public:
   Bundle(std::vector<object::fiber::Fiber> fibers,
          std::vector<layer::Property> properties);
   ~Bundle() = default;

   // getter
   const Fiber &fiber(size_t i) const { return fibers_[i]; }
   const std::vector<Fiber> &fibers() const { return fibers_; }
   const aabb::AABB<float, 3> &voi() const { return voi_; }
   layer::Properties layers() const { return layers_; };
   size_t size() const { return fibers_.size(); }
   bool empty() const { return fibers_.empty(); }

   // layer informations
   size_t layer_size() const { return layers_.size(); }
   float layer_scale(size_t i) const { return layers_.scale()[i]; }
   float layer_dn(size_t i) const { return layers_.dn()[i]; }
   float layer_mu(size_t i) const { return layers_.mu()[i]; }
   layer::Orientation layer_orientation(size_t i) const {
      return layers_.orientation()[i];
   }
   const std::vector<float> &layers_scale() const { return layers_.scale(); }
   const std::vector<float> &layers_dn() const { return layers_.dn(); }
   const std::vector<float> &layers_mu() const { return layers_.mu(); }
   const std::vector<layer::Orientation> &layers_orientation() const {
      return layers_.orientation();
   }
   const std::vector<float> &layers_scale_sqr() const {
      return layers_.scale_sqr();
   }

   // manipulator
   void Resize(const float f);
   void ResizePoints(const float f);
   void ResizeRadii(const float f);
   void Rotate(const vm::Mat3x3<float> &rot_mat);
   void RotateAroundPoint(const vm::Mat3x3<float> &rot_mat,
                          const vm::Vec3<float> &point);
   void Translate(const vm::Vec3<float> &translation);

 private:
   std::vector<Fiber> fibers_{};
   layer::Properties layers_{};
   aabb::AABB<float, 3> voi_{};

   void CalculateVoi();
};
} // namespace fiber

#endif // SIMULATION_FIBER_BUNDLE_HPP_
