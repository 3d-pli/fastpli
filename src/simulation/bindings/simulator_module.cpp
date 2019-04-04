#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../simulator.hpp"
#include "objects/np_array_container.hpp"
#include "objects/np_array_helper.hpp"

namespace py = pybind11;

PYBIND11_MODULE(__simulation, m) {
   m.doc() = "Simulation of 3D-PLI images";

   py::class_<PliSimulator>(m, "_Simulator")
       .def(py::init())
       .def("set_pli_setup", &PliSimulator::SetPliSetup)
       .def("set_tissue_properties", &PliSimulator::SetTissueProperties)
       .def("run_simulation",
            [](PliSimulator &self, std::array<long long, 3> dim,
               py::array_t<int, py::array::c_style> label_array,
               py::array_t<float, py::array::c_style> vector_array,
               double theta, double phi, double ps, bool do_nn = true) {
               auto label_container = object::NpArray2Container(label_array);
               auto vector_container = object::NpArray2Container(vector_array);

               auto image = new std::vector<float>(
                   self.RunSimulation(dim, label_container, vector_container,
                                      theta, phi, ps, do_nn));

               std::vector<size_t> dim_image =
                   vm::cast<size_t>(self.GetImageDim());

               return object::Vec2NpArray(image, dim_image);
            },
            py::arg("dim"), py::arg("label_field"), py::arg("vector_field"),
            py::arg("theta") = 0, py::arg("phi") = 0, py::arg("step_size") = 1,
            py::arg("do_nn_intp") = true);

   py::class_<PliSimulator::PhyProp>(m, "_PhyProp")
       .def(py::init<double, double>())
       .def_readwrite("dn", &PliSimulator::PhyProp::dn)
       .def_readwrite("mu", &PliSimulator::PhyProp::mu);

   py::class_<PliSimulator::Setup>(m, "_Setup")
       .def(py::init())
       .def_readwrite("light_intensity", &PliSimulator::Setup::light_intensity)
       .def_readwrite("pixel_size", &PliSimulator::Setup::pixel_size)
       .def_readwrite("wavelength", &PliSimulator::Setup::wavelength)
       .def_readwrite("untilt_sensor", &PliSimulator::Setup::untilt_sensor)
       .def_readwrite("filter_rotations",
                      &PliSimulator::Setup::filter_rotations);
}
