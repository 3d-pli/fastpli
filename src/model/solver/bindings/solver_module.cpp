#include <tuple>

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../world.hpp"
#include "objects/fiber.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_solver_cpp, m) {
   m.doc() = "Fiber Volume Colliding Solver";

   py::class_<World>(m, "_SolverCPP")
       .def(py::init())
       //  .def_property("fiber_bundles", &World::get_fibers,
       //  &World::set_fibers)
       .def("get_fiber_bundles", &World::get_fibers)
       .def("set_fiber_bundles", &World::set_fibers)
       .def("set_parameters",
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
       .def("get_parameters",
            [](World &self) {
               auto p = self.get_parameter();
               return std::make_tuple(p.drag, p.obj_min_radius,
                                      p.obj_mean_length);
            })
       .def("step", (bool (World::*)(void)) & World::Step)
       .def_property_readonly("num_objects", &World::NumObj);
}