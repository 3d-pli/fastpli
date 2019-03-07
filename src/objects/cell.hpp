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
   float CalcRadius(size_t idx, float t) const = delete;
};

using CellPopulation = std::vector<Cell>;
using CellPopulations = std::vector<CellPopulation>;

} // namespace object

#endif // CELL_HPP_
