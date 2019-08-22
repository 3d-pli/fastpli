#include "cell.hpp"

namespace object {
Cell::Cell(const std::vector<double> &points,
           const std::vector<double> &radii) {

   if (points.size() != radii.size() * 3)
      throw std::invalid_argument("points and radii aren't same size");

   points_.resize(radii.size());
   points_.shrink_to_fit();
   voi_ = aabb::AABB<double, 3>{};
   radii_ = radii;

   for (size_t i = 0; i < radii_.size(); i++)
      points_[i] = vm::Vec3<double>(points[3 * i + 0], points[3 * i + 1],
                                    points[3 * i + 2]);

   // calc voi
   if (points_.empty())
      return;
   else if (points_.size() == 1) {
      voi_ = aabb::AABB<double, 3>(points_[0], points_[0]);
      return;
   }
   CalculateVoi();
}

Cell::Cell(const std::vector<vm::Vec3<double>> &points,
           const std::vector<double> &radii) {

   if (points.size() != radii.size())
      throw std::invalid_argument("points and radii aren't same size");

   points_ = points;
   radii_ = radii;

   // calc voi
   if (points_.empty())
      return;
   else if (points_.size() == 1) {
      voi_ = aabb::AABB<double, 3>(points_[0], points_[0]);
      return;
   }
   CalculateVoi();
}
} // namespace object
