#ifndef SIMULATION_CELL_POPULATION_HPP_
#define SIMULATION_CELL_POPULATION_HPP_

#include <utility>
#include <vector>

#include "include/aabb.hpp"
#include "include/vemath.hpp"
#include "objects/cell.hpp"

namespace cell {

class Property {
 public:
   Property() = default;
   Property(double s, double m) : scale_(s), mu_(m) { scale_sqr_ = s * s; }

   double scale() const { return scale_; }
   double scale_sqr() const { return scale_sqr_; }
   double mu() const { return mu_; }

 private:
   double scale_{0};
   double scale_sqr_{0};
   double mu_{0};
};

class Population {
 public:
   Population(std::vector<std::vector<double>> cells, Property property);
   Population(std::vector<object::Cell> cells, Property property);
   ~Population() = default;

   // getter
   const object::Cell &cell(size_t i) const { return cells_[i]; }
   const std::vector<object::Cell> &cells() const { return cells_; }
   const aabb::AABB<double, 3> &voi() const { return voi_; }
   size_t size() const { return cells_.size(); }
   bool empty() const { return cells_.empty(); }

   // properties
   double scale() const { return property_.scale(); }
   double scale_sqr() const { return property_.scale_sqr(); }
   double mu() const { return property_.mu(); }

   // manipulator
   void Resize(const double f);
   void ResizePoints(const double f);
   void ResizeRadii(const double f);
   void Rotate(const vm::Mat3x3<double> &rot_mat);
   void RotateAroundPoint(const vm::Mat3x3<double> &rot_mat,
                          const vm::Vec3<double> &point);
   void Translate(const vm::Vec3<double> &translation);

 private:
   std::vector<object::Cell> cells_{};
   Property property_{};
   aabb::AABB<double, 3> voi_{};

   void CalculateVoi();
};
} // namespace cell

#endif // SIMULATION_CELL_POPULATION_HPP_
