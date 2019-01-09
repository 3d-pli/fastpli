#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../pli_generator.hpp"
#include "objects/vector_container.hpp"

namespace py = pybind11;

PYBIND11_MODULE(generation, m) {
   m.doc() = "Generation of tissue for SimPLI";

   py::class_<PliGenerator>(m, "Generator")
       .def(py::init())
       .def("set_volume", &PliGenerator::SetVolumeWithArrays, py::arg("dim"),
            py::arg("pixel_size"),
            py::arg("flip_direction") =
                std::array<bool, 3>{{false, false, false}})
       .def("set_fiber_bundles",
            [](PliGenerator &self, std::vector<std::vector<fiber::Data>> fbs,
               std::vector<std::vector<fiber::layer::Property>> prs) {
               if (fbs.size() != prs.size())
                  throw py::value_error("fbs and prs not the same size");

               std::vector<fiber::Bundle> fiber_bundles;
               for (size_t i = 0; i < fbs.size(); i++)
                  fiber_bundles.emplace_back(fiber::Bundle(fbs[i], prs[i]));

               self.SetFiberBundles(fiber_bundles);
            })
       .def("run_generation", &PliGenerator::RunTissueGeneration,
            py::arg("only_label") = false, py::arg("debug") = false);

   // TODO: to objects?
   py::enum_<fiber::layer::Orientation>(m, "Orientation")
       .value("background", fiber::layer::Orientation::background)
       .value("parallel", fiber::layer::Orientation::parallel)
       .value("radial", fiber::layer::Orientation::radial);

   py::class_<fiber::layer::Property>(m, "LayerProperty")
       .def(py::init<float, float, float, char>())
       .def(py::init<float, float, float, fiber::layer::Orientation>())
       .def_readwrite("scale", &fiber::layer::Property::scale)
       .def_readwrite("dn", &fiber::layer::Property::dn)
       .def_readwrite("mu", &fiber::layer::Property::mu)
       .def_readwrite("layer_orientation",
                      &fiber::layer::Property::orientation);
}