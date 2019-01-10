#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../simulator.hpp"
#include "objects/vector_container.hpp"

namespace py = pybind11;

PYBIND11_MODULE(simulation, m) {
   m.doc() = "Simulation of 3D-PLI images";

   py::class_<PliSimulator>(m, "Simulator")
       .def(py::init())
       .def("set_pli_setup", &PliSimulator::SetPliSetup)
       .def("set_tissue",
            (void (PliSimulator::*)(
                data::VectorContainer<int>, data::VectorContainer<float>,
                const std::array<int, 3> &,
                const std::vector<PliSimulator::TissueProperty> &,
                const double)) &
                PliSimulator::SetTissue)
       .def("run_simulation",
            [](PliSimulator &self, double theta, double phi, double ps,
               bool do_nn = true) {
               auto result = new std::vector<float>(
                   self.RunSimulation(theta, phi, ps, do_nn));
               auto capsule = py::capsule(result, [](void *result) {
                  delete reinterpret_cast<std::vector<float> *>(result);
               });
               return py::array(result->size(), result->data(), capsule);
            },
            py::arg("theta") = 0, py::arg("phi") = 0, py::arg("step_size") = 1,
            py::arg("do_nn_intp") = true);

   py::class_<PliSimulator::TissueProperty>(m, "TissueProperty")
       .def(py::init())
       .def_readwrite("dn", &PliSimulator::TissueProperty::dn)
       .def_readwrite("mu", &PliSimulator::TissueProperty::mu);

   py::class_<PliSimulator::Setup>(m, "Setup")
       .def(py::init())
       .def_readwrite("light_intensity", &PliSimulator::Setup::light_intensity)
       .def_readwrite("resolution", &PliSimulator::Setup::resolution)
       .def_readwrite("wavelength", &PliSimulator::Setup::wavelength)
       .def_readwrite("untilt_sensor", &PliSimulator::Setup::untilt_sensor)
       .def_readwrite("filter_rotations",
                      &PliSimulator::Setup::filter_rotations);
}
