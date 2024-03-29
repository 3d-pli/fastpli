#ifndef SRC_SIMULATION_FIBER_BUNDLE_HPP_
#define SRC_SIMULATION_FIBER_BUNDLE_HPP_

#include <utility>
#include <vector>

#include "include/aabb.hpp"
#include "include/vemath.hpp"
#include "objects/fiber.hpp"

namespace fiber {

namespace layer {
enum class Orientation { none, parallel, radial };

struct Property {
   Property(double s, double n, double m, char o);
   Property(double s, double n, double m, Orientation o);

   double scale{};
   double dn{};
   double mu{};
   Orientation orientation{};
};

class Properties {
 public:
   const std::vector<double> &scale_squ() const { return scale_squ_; }
   const std::vector<double> &dn() const { return dn_; }
   const std::vector<double> &mu() const { return mu_; }
   const std::vector<Orientation> &orientation() const { return orientation_; }

   size_t size() const;
   void clear();
   void reserve(size_t i);
   void resize(size_t i);
   void push_back(Property p);
   void push_back(double s, double n, double m, char o) {
      push_back(Property(s, n, m, o));
   }
   void push_back(double s, double n, double m, Orientation o) {
      push_back(Property(s, n, m, o));
   }

 private:
   std::vector<double> scale_squ_;
   std::vector<double> dn_;
   std::vector<double> mu_;
   std::vector<Orientation> orientation_;
};

} // namespace layer

class Bundle {
 public:
   Bundle(std::vector<object::Fiber> fibers,
          std::vector<layer::Property> properties);
   Bundle(std::vector<std::vector<double>> fibers,
          std::vector<layer::Property> properties);

   // getter
   const object::Fiber &fiber(size_t i) const { return fibers_[i]; }
   const std::vector<object::Fiber> &fibers() const { return fibers_; }
   const aabb::AABB<double, 3> &aabb() const { return aabb_; }
   layer::Properties layers() const { return layers_; }
   size_t size() const { return fibers_.size(); }
   bool empty() const { return fibers_.empty(); }

   // layer informations
   size_t layer_size() const { return layers_.size(); }
   double layer_dn(size_t i) const { return layers_.dn()[i]; }
   double layer_mu(size_t i) const { return layers_.mu()[i]; }
   layer::Orientation layer_orientation(size_t i) const {
      return layers_.orientation()[i];
   }
   const std::vector<double> &layers_dn() const { return layers_.dn(); }
   const std::vector<double> &layers_mu() const { return layers_.mu(); }
   const std::vector<layer::Orientation> &layers_orientation() const {
      return layers_.orientation();
   }
   const std::vector<double> &layers_scale_squ() const {
      return layers_.scale_squ();
   }

   // manipulator
   void Resize(const double f);
   void ResizePoints(const double f);
   void ResizeRadii(const double f);
   void Rotate(const vm::Mat3x3<double> &rot_mat);
   void RotateAroundPoint(const vm::Mat3x3<double> &rot_mat,
                          const vm::Vec3<double> &point);
   void Translate(const vm::Vec3<double> &translation);

 private:
   object::FiberBundle fibers_{};
   layer::Properties layers_{};
   aabb::AABB<double, 3> aabb_{};

   void CalculateVoi();
};
} // namespace fiber

#endif // SRC_SIMULATION_FIBER_BUNDLE_HPP_
