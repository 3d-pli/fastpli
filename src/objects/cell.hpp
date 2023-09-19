#ifndef SRC_OBJECTS_CELL_HPP_
#define SRC_OBJECTS_CELL_HPP_

#include <cassert>
#include <vector>

#include "include/aabb.hpp"
#include "include/vemath.hpp"
#include "objects/fiber.hpp"

namespace object {

class Cell : public object::Fiber {
 public:
   Cell(const std::vector<vm::Vec3<double>> &points,
        const std::vector<double> &radii);
   Cell(const std::vector<double> &points, const std::vector<double> &radii);

   // defaults
   Cell() = default;
   Cell(Cell &&) = default;
   Cell(const Cell &) = default;
   Cell &operator=(Cell &&) = default;
   Cell &operator=(const Cell &) = default;

   ~Cell() = default;

   double CalcRadius(size_t idx, double t) const = delete;
};

using CellPopulation = std::vector<Cell>;
using CellPopulations = std::vector<CellPopulation>;

} // namespace object

#endif // SRC_OBJECTS_CELL_HPP_
