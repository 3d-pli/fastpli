#ifndef FIBER_MODEL_HPP_
#define FIBER_MODEL_HPP_

#include <utility>
#include <vector>

#include "include/aabb.hpp"
#include "include/vemath.hpp"
#include "objects/fiber.hpp"

enum class Orientation { background, parallel, radial };

class FiberData {
 public:
   FiberData(const std::vector<vm::Vec3<float>> &points,
             const std::vector<float> &radii);
   FiberData(const std::vector<float> &points,
             const std::vector<float> &radii);
   FiberData(const data::Fiber &fiber_data)
       : FiberData(fiber_data.points(), fiber_data.radii()){};

   // defaults
   FiberData() = default;
   FiberData(FiberData &&) = default;
   FiberData(const FiberData &) = default;
   FiberData &operator=(FiberData &&) = default;
   FiberData &operator=(const FiberData &) = default;

   ~FiberData() = default;

   // getter
   // TODO: return of arrray should be only in python module
   const std::vector<vm::Vec3<double>> &pos() const { return pos_; }
   const std::vector<double> &radii() const { return radii_; }
   const aabb::AABB<double, 3> &voi() const { return voi_; }

   size_t size() const { return pos_.size(); };
   double CalcRadius(size_t idx, double t) const;

   // manipolator
   void RotateFiber(const std::array<double, 9> &rot_mat) {
      RotateFiber(vm::Mat3x3<double>(rot_mat));
   };
   void RotateFiber(const vm::Mat3x3<double> &rot_mat);

   void RotateFiberAroundPoint(const std::array<double, 9> &rot_mat,
                               const std::array<double, 3> &point) {
      RotateFiberAroundPoint(vm::Mat3x3<double>(rot_mat),
                             vm::Vec3<double>(point));
   };

   void RotateFiberAroundPoint(const vm::Mat3x3<double> &rot_mat,
                               const vm::Vec3<double> point);

   void TranslateFiber(const std::array<double, 3> &translation) {
      TranslateFiber(vm::Vec3<double>(translation));
   };

   void TranslateFiber(const vm::Vec3<double> &translation);

   void ResizeFiberPos(const double f);
   void ResizeFiberRadii(const double f);
   void ResizeFiber(const double f);

 private:
   std::vector<vm::Vec3<double>> pos_;
   std::vector<double> radii_;
   aabb::AABB<double, 3> voi_;

   void CalculateVoi();
};

struct LayerProperty {
   LayerProperty() {}
   LayerProperty(double s, double n, double m, ushort o)
       : scale(s), dn(n), mu(m) {
      if (o == 0)
         orientation = Orientation::background;
      else if (o == 1)
         orientation = Orientation::parallel;
      else if (o == 2)
         orientation = Orientation::radial;
      else
         throw std::invalid_argument("Orientation must be 0, 1, or 2");
   }
   LayerProperty(double s, double n, double m, Orientation o)
       : scale(s), dn(n), mu(m), orientation(o) {}
   LayerProperty(std::array<double, 4> data)
       : scale(data[0]), dn(data[1]), mu(data[2]) {
      if (data[3] == 0)
         orientation = Orientation::background;
      else if (data[3] == 1)
         orientation = Orientation::parallel;
      else if (data[3] == 2)
         orientation = Orientation::radial;
      else
         throw std::invalid_argument("Orientation must be 0, 1, or 2");
   }

   double scale{};
   double dn{};
   double mu{};
   Orientation orientation{Orientation::background};
};

class FiberBundle {
 public:
   FiberBundle() = default;
   ~FiberBundle() = default;

   // setter
   void push_fiber(FiberData &fiber);
   void SetFiberBundleProperties(const std::vector<LayerProperty> &prop);

   // manipolator
   void RotateFiberBundle(const std::array<double, 9> &rot_mat);
   void RotateFiberBundleAroundPoint(const std::array<double, 9> &rot_mat,
                                     std::array<double, 3> point);
   void TranslateFiberBundle(const std::array<double, 3> &translation);

   void ResizeFiberBundlePos(const double f);
   void ResizeFiberBundleRadii(const double f);
   void ResizeFiberBundle(const double f);

   // getter
   const FiberData &fiber(size_t i) const { return fibers_[i]; }
   const std::vector<FiberData> &fibers() const { return fibers_; }
   size_t size() const { return fibers_.size(); }
   const std::vector<double> &layer_radii() const { return layer_scale_; }
   const std::vector<double> &layer_scale_sqr() const {
      return layer_scale_sqr_;
   }
   const std::vector<Orientation> &layer_orientation() const {
      return layer_orientation_;
   }

   size_t layer_size() const { return layer_dn_.size(); }
   double layer_dn(size_t i) const { return layer_dn_[i]; }
   double layer_mu(size_t i) const { return layer_mu_[i]; }
   const std::vector<double> &layer_dn() const { return layer_dn_; }
   const std::vector<double> &layer_mu() const { return layer_mu_; }
   const aabb::AABB<double, 3> &voi() const { return voi_; }

 private:
   std::vector<FiberData> fibers_;
   aabb::AABB<double, 3> voi_{};

   std::vector<double> layer_scale_;
   std::vector<double> layer_scale_sqr_;
   std::vector<double> layer_dn_;
   std::vector<double> layer_mu_;
   std::vector<Orientation> layer_orientation_;
};

#endif // FIBER_MODEL_HPP_