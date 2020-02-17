#include <iostream>
#include <tuple>

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../world.hpp"
#include "objects/fiber.hpp"
#include "objects/np_array_container.hpp"
#include "objects/np_array_helper.hpp"

namespace py = pybind11;

PYBIND11_MODULE(__solver, m) {
   m.doc() = "Fiber Volume Colliding Solver";

   py::class_<World>(m, "_Solver")
       .def(py::init())
       .def("_set_fiber_bundles",
            [](World &self,
               std::vector<std::vector<py::array_t<double, py::array::c_style>>>
                   fbs) {
               std::vector<std::vector<std::vector<double>>> fiber_bundles;
               for (auto i = 0u; i < fbs.size(); i++) {
                  auto fiber_bundle = std::vector<std::vector<double>>();
                  for (auto &f : fbs[i]) {
                     fiber_bundle.push_back(object::NpArray2Vector<double>(f));
                  }
                  fiber_bundles.push_back(fiber_bundle);
               }
               self.set_fibers_vector(fiber_bundles);
            })
       .def("_get_fiber_bundles",
            [](World &self) {
               auto const fbs = self.get_fibers_vector();
               std::vector<size_t> dim{0, 4};
               std::vector<std::vector<py::array>> fiber_bundles;
               for (size_t i = 0; i < fbs.size(); i++) {
                  if (fbs[i].size() > 0)
                     fiber_bundles.push_back(std::vector<py::array>());

                  for (size_t j = 0; j < fbs[i].size(); j++) {
                     dim[0] = fbs[i][j].size() / 4;
                     fiber_bundles[i].push_back(object::Vec2NpArray(
                         new std::vector<double>(fbs[i][j]), dim));
                  }
               }

               return fiber_bundles;
            })
       .def("_set_omp_num_threads", &World::set_omp_num_threads)
       .def(
           "_set_parameters",
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
       .def("_get_parameters",
            [](World &self) {
               auto p = self.get_parameter();
               return std::make_tuple(p.drag, p.obj_min_radius,
                                      p.obj_mean_length);
            })
       .def("_set_col_voi",
            [](World &self, std::array<double, 3> min,
               std::array<double, 3> max) {
               self.set_colliding_voi(aabb::AABB<double, 3>(min, max));
            })
       .def("step", (bool (World::*)(void)) & World::Step)
       .def("apply_boundary_conditions", &World::ApplyBoundaryConditions,
            py::arg("steps") = 1)
       .def_property_readonly("num_obj", &World::num_obj)
       .def_property_readonly("num_col_obj", &World::num_col_obj)
       .def_property_readonly("overlap", &World::overlap)
       .def_property_readonly("max_speed", &World::max_speed)
       .def("draw_scene", &World::DrawScene, py::arg("rot_x") = 30,
            py::arg("rot_y") = 30, py::arg("rot_z") = 0,
            py::arg("only_col") = false)
       .def("save_ppm", &World::SavePPM, py::arg("file"));
}
