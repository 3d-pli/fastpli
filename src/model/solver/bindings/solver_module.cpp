#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../world.hpp"
// #include "fiber.hpp"
#include "include/vemath.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_model, m) {
   m.doc() = "Volume Colliding Solver World";

   py::class_<World>(m, "Solver")
       .def(py::init())
       .def("get_fibers", &World::get_fibers)
       .def("set_fibers", &World::set_fibers)
       .def("set_parameter",
            [](World &self, float drag, float obj_min_radius,
               float obj_mean_length) {
               World::WorldParameter wp;
               wp.drag = drag;
               wp.obj_min_radius = obj_min_radius;
               wp.obj_mean_length = obj_mean_length;
               self.set_parameter(wp);
            },
            py::arg("drag") = 0, py::arg("obj_min_radius") = 10,
            py::arg("obj_mean_length") = 10)
       .def("step", (bool (World::*)(void)) & World::Step)
       .def_property_readonly("num_objects", &World::NumObj);
}