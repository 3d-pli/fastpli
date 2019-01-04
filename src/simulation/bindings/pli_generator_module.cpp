#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../core/pli_generator.hpp"
#include "pli_module.hpp"

namespace py = pybind11;

PYBIND11_MODULE(generation, m) {
   m.doc() = "Generation of tissue for SimPLI";

   py::class_<PliGenerator>(m, "Generator")
       .def(py::init())
       .def("set_volume", &PliGenerator::SetVolumeWithArrays, py::arg("dim"),
            py::arg("pixel_size"),
            py::arg("flip_direction") =
                std::array<bool, 3>{{false, false, false}})
       .def("set_fiber_bundle", &PliGenerator::SetFiberBundles)
       .def("rotate_fiber_bundle", &PliGenerator::RotateFiberBundles)
       .def("translate_fiber_bundle", &PliGenerator::TranslateFiberBundles)
       .def("rotate_fiber_around_fiber_bundle",
            &PliGenerator::RotateFiberBundlesAroundPoint)
       .def("run_generation", &PliGenerator::RunTissueGeneration,
            py::arg("only_label") = false, py::arg("debug") = false);

   py::enum_<Orientation>(m, "Orientation")
       .value("background", Orientation::background)
       .value("parallel", Orientation::parallel)
       .value("radial", Orientation::radial);
   //  .export_values();

   py::class_<FiberData>(m, "FiberData")
       .def(
           py::init<const std::vector<double> &, const std::vector<double> &>())
       // TODO: numpy 2d array
       .def("pos", &FiberData::pos)
       .def("radii", &FiberData::radii)
       .def("size", &FiberData::size)
       .def("calc_radius", &FiberData::CalcRadius)
       .def("rotate_fiber",
            (void (FiberData::*)(const std::array<double, 9> &)) &
                FiberData::RotateFiber)
       .def("rotate_fiber_around_point",
            (void (FiberData::*)(const std::array<double, 9> &,
                                 const std::array<double, 3> &)) &
                FiberData::RotateFiberAroundPoint)
       .def("translate_fiber",
            (void (FiberData::*)(const std::array<double, 3> &)) &
                FiberData::TranslateFiber)
       .def("resize_fiber_pos", &FiberData::ResizeFiberPos)
       .def("resize_fiber_radii", &FiberData::ResizeFiberRadii)
       .def("resize_fiber", &FiberData::ResizeFiber);

   py::class_<FiberBundle>(m, "FiberBundle")
       .def(py::init())
       .def("push_fiber", &FiberBundle::push_fiber)
       .def("set_fiber_bundle_properties",
            &FiberBundle::SetFiberBundleProperties)
       .def("rotate_fiber_bundle", &FiberBundle::RotateFiberBundle)
       .def("rotate_fiber_bundle_around_point",
            &FiberBundle::RotateFiberBundleAroundPoint)
       .def("translate_fiber_bundle", &FiberBundle::TranslateFiberBundle)
       .def("resize_fiber_bundle_pos", &FiberBundle::ResizeFiberBundlePos)
       .def("resize_fiber_bundle_radii", &FiberBundle::ResizeFiberBundleRadii)
       .def("resize_fiber_bundle", &FiberBundle::ResizeFiberBundle);

   py::class_<LayerProperty>(m, "LayerProperty")
       .def(py::init<double, double, double, int>())
       .def(py::init<double, double, double, Orientation>())
       .def(py::init<std::array<double, 4>>())
       .def_readwrite("scale", &LayerProperty::scale)
       .def_readwrite("dn", &LayerProperty::dn)
       .def_readwrite("mu", &LayerProperty::mu)
       .def_readwrite("layer_orientation", &LayerProperty::orientation);
}
