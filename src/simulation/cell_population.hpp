#ifndef SIMULATION_CELL_POPULATION_HPP_
#define SIMULATION_CELL_POPULATION_HPP_

#include <utility>
#include <vector>

#include "include/aabb.hpp"
#include "include/vemath.hpp"
#include "objects/cell.hpp"

namespace cell {

struct Property {
   Property() = default;
   Property(float m) : mu(m) {}

   // could be more
   float mu{};
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
   float mu() const { return property_.mu; }

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
