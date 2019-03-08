#ifndef CELL_HPP_
#define CELL_HPP_

#include <cassert>
#include <vector>

#include "include/aabb.hpp"
#include "include/vemath.hpp"
#include "objects/fiber.hpp"

namespace object {

class Cell : public object::Fiber {
 public:
   Cell(const std::vector<vm::Vec3<float>> &points,
        const std::vector<float> &radii);
   Cell(const std::vector<float> &points, const std::vector<float> &radii);

   // defaults
   Cell() = default;
   Cell(Cell &&) = default;
   Cell(const Cell &) = default;
   Cell &operator=(Cell &&) = default;
   Cell &operator=(const Cell &) = default;

   ~Cell() = default;

   float CalcRadius(size_t idx, float t) const = delete;
};

using CellPopulation = std::vector<Cell>;
using CellPopulations = std::vector<CellPopulation>;

} // namespace object

#endif // CELL_HPP_
