
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../fiber.hpp"

namespace py = pybind11;

void fiber_class(py::module &);
void cell_class(py::module &);

PYBIND11_MODULE(__objects, m) {
   fiber_class(m);
   cell_class(m);
}
