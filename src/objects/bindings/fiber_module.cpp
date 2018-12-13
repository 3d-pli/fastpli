#include "../core/fiber.hpp"

#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(_fiber_cpp, m) { m.doc() = "fiber data as c++ class"; }
