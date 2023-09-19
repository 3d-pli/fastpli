#ifndef SRC_SIMULATION_CELL_POPULATION_HPP_
#define SRC_SIMULATION_CELL_POPULATION_HPP_

#include <utility>
#include <vector>

#include "include/aabb.hpp"
#include "include/vemath.hpp"
#include "objects/cell.hpp"

namespace cell {

class Property {
 public:
   Property() = default;
   Property(double s, double m) : scale_(s), mu_(m) { scale_squ_ = s * s; }

   double scale() const { return scale_; }
   double scale_squ() const { return scale_squ_; }
   double mu() const { return mu_; }

 private:
   double scale_{0};
   double scale_squ_{0};
   double mu_{0};
};

class Population {
 public:
   Population(std::vector<std::vector<double>> cells, Property property);
   Population(std::vector<object::Cell> cells, Property property);

   // getter
   const object::Cell &cell(size_t i) const { return cells_[i]; }
   const std::vector<object::Cell> &cells() const { return cells_; }
   const aabb::AABB<double, 3> &aabb() const { return aabb_; }
   size_t size() const { return cells_.size(); }
   bool empty() const { return cells_.empty(); }

   // properties
   double scale() const { return property_.scale(); }
   double scale_squ() const { return property_.scale_squ(); }
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
   aabb::AABB<double, 3> aabb_{};

   void CalculateVoi();
};
} // namespace cell

#endif // SRC_SIMULATION_CELL_POPULATION_HPP_
