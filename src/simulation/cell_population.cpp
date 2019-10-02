#include "cell_population.hpp"

#include <cassert>
#include <iostream>
#include <vector>

#include "include/aabb.hpp"
#include "include/vemath.hpp"
#include "objects/cell.hpp"

namespace cell {

Population::Population(std::vector<std::vector<double>> cells,
                       Property property) {

   for (auto const c : cells) {
      std::vector<vm::Vec3<double>> points;
      std::vector<double> radii;

      assert(c.size() % 4 == 0);

      points.reserve(c.size() / 4);
      radii.reserve(c.size() / 4);
      for (auto i = 0u; i < c.size(); i = i + 4) {
         points.push_back(vm::Vec3<double>(c[i + 0], c[i + 1], c[i + 2]));
         radii.push_back(c[i + 3]);
      }

      cells_.emplace_back(object::Cell(points, radii));
   }
   property_ = property;
   CalculateVoi();
}

Population::Population(std::vector<object::Cell> cells, Property property) {
   cells_.swap(cells);
   property_ = property;
   CalculateVoi();
}

void Population::Resize(const double f) {
   for (auto &cell : cells_)
      cell.Resize(f);
   CalculateVoi();
}
void Population::ResizePoints(const double f) {
   for (auto &cell : cells_)
      cell.ResizePoints(f);
   CalculateVoi();
}

void Population::ResizeRadii(const double f) {
   for (auto &cell : cells_)
      cell.ResizeRadii(f);
   CalculateVoi();
}

void Population::Rotate(const vm::Mat3x3<double> &rot_mat) {
   for (auto &cell : cells_)
      cell.Rotate(rot_mat);
   CalculateVoi();
}

void Population::RotateAroundPoint(const vm::Mat3x3<double> &rot_mat,
                                   const vm::Vec3<double> &point) {
   for (auto &cell : cells_)
      cell.RotateAroundPoint(rot_mat, point);
   CalculateVoi();
}

void Population::Translate(const vm::Vec3<double> &translation) {
   for (auto &cell : cells_)
      cell.Translate(translation);
   CalculateVoi();
}

void Population::CalculateVoi() {
   voi_ = aabb::AABB<double, 3>();
   if (!cells_.empty())
      voi_ = cells_.front().voi();

   for (auto &cell : cells_)
      voi_.Unite(cell.voi());
}
} // namespace cell
