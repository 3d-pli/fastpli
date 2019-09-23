#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <iostream>
#include <mpi.h>

#include "../setup.hpp"
#include "../simulator.hpp"
#include "objects/np_array_container.hpp"
#include "objects/np_array_helper.hpp"

namespace py = pybind11;

PYBIND11_MODULE(__simulation, m) {
   m.doc() = "Simulation of 3D-PLI images";

   py::class_<PliSimulator>(m, "_Simulator")
       .def(py::init())
       .def(
           "set_pli_setup",
           [](PliSimulator &self, double step_size, double light_intensity,
              double voxel_size, double wavelength, bool interpolate,
              bool untilt_sensor, bool flip_z,
              std::vector<double> filter_rotations) {
              self.SetSetup(setup::Setup(step_size, light_intensity, voxel_size,
                                         wavelength, interpolate, untilt_sensor,
                                         flip_z, filter_rotations));
           },
           py::arg("step_size"), py::arg("light_intensity"),
           py::arg("voxel_size"), py::arg("wavelength"), py::arg("interpolate"),
           py::arg("untilt_sensor"), py::arg("flip_z"),
           py::arg("filter_rotations"))
       .def("set_mpi_comm",
            [](PliSimulator &self, long long comm_address) {
               MPI_Comm comm = *static_cast<MPI_Comm *>(
                   reinterpret_cast<void *>(comm_address));
               self.SetMPIComm(comm);
            })
       .def(
           "run_simulation",
           [](PliSimulator &self, std::array<long long, 3> dim,
              py::array_t<int, py::array::c_style> label_array,
              py::array_t<float, py::array::c_style> vector_array,
              py::array_t<double, py::array::c_style> prop_array, double theta,
              double phi) {
              auto label_container = object::NpArray2Container(label_array);
              auto vector_container = object::NpArray2Container(vector_array);
              auto properties =
                  setup::PhyProps(object::NpArray2Vector<double>(prop_array));

              auto image = new std::vector<double>(
                  self.RunSimulation(dim, label_container, vector_container,
                                     properties, theta, phi));

              std::vector<size_t> dim_image =
                  vm::cast<size_t>(self.GetImageDim());

              return object::Vec2NpArray(image, dim_image);
           },
           py::arg("dim"), py::arg("label_field"), py::arg("vector_field"),
           py::arg("properties"), py::arg("theta") = 0, py::arg("phi") = 0)
       .def("set_omp_num_threads", &PliSimulator::set_omp_num_threads);
}
