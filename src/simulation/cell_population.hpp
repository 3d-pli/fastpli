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
   Property(float s, float m) : scale_(s), mu_(m) { scale_sqr_ = s * s; }

   float scale() const { return scale_; }
   float scale_sqr() const { return scale_sqr_; }
   float mu() const { return mu_; }

 private:
   float scale_{0};
   float scale_sqr_{0};
   float mu_{0};
};

class Population {
 public:
   Population(std::vector<object::Cell> cells, Property property);
   ~Population() = default;

   // getter
   const object::Cell &cell(size_t i) const { return cells_[i]; }
   const std::vector<object::Cell> &cells() const { return cells_; }
   const aabb::AABB<float, 3> &voi() const { return voi_; }
   size_t size() const { return cells_.size(); }
   bool empty() const { return cells_.empty(); }

   // properties
   float scale() const { return property_.scale(); }
   float scale_sqr() const { return property_.scale_sqr(); }
   float mu() const { return property_.mu(); }

   // manipulator
   void Resize(const float f);
   void ResizePoints(const float f);
   void ResizeRadii(const float f);
   void Rotate(const vm::Mat3x3<float> &rot_mat);
   void RotateAroundPoint(const vm::Mat3x3<float> &rot_mat,
                          const vm::Vec3<float> &point);
   void Translate(const vm::Vec3<float> &translation);

 private:
   std::vector<object::Cell> cells_{};
   Property property_{};
   aabb::AABB<float, 3> voi_{};

   void CalculateVoi();
};
} // namespace cell

#endif // SIMULATION_CELL_POPULATION_HPP_
