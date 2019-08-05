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
       .def("set_pli_setup",
            [](PliSimulator &self, double light_intensity, double voxel_size,
               double wavelength, std::vector<double> filter_rotations,
               bool untilt_sensor) {
               setup::Setup setup;
               setup.light_intensity = light_intensity;
               setup.wavelength = wavelength;
               setup.voxel_size = voxel_size;
               setup.untilt_sensor = untilt_sensor;
               setup.filter_rotations = filter_rotations;
               self.SetPliSetup(setup);
            },
            py::arg("light_intensity"), py::arg("wavelength"),
            py::arg("voxel_size"), py::arg("untilt_sensor"),
            py::arg("filter_rotations"))
       .def("set_mpi_comm",
            [](PliSimulator &self, long long comm_address) {
               MPI_Comm comm = *static_cast<MPI_Comm *>(
                   reinterpret_cast<void *>(comm_address));
               self.SetMPIComm(comm);
            })
       .def("run_simulation",
            [](PliSimulator &self, std::array<long long, 3> dim,
               py::array_t<int, py::array::c_style> label_array,
               py::array_t<float, py::array::c_style> vector_array,
               py::array_t<double, py::array::c_style> prop_array, double theta,
               double phi, double ps, bool do_nn = true) {
               auto label_container = object::NpArray2Container(label_array);
               auto vector_container = object::NpArray2Container(vector_array);
               auto properties =
                   setup::PhyProps(object::NpArray2Vector<double>(prop_array));

               auto image = new std::vector<float>(
                   self.RunSimulation(dim, label_container, vector_container,
                                      properties, theta, phi, ps, do_nn));

               std::vector<size_t> dim_image =
                   vm::cast<size_t>(self.GetImageDim());

               return object::Vec2NpArray(image, dim_image);
            },
            py::arg("dim"), py::arg("label_field"), py::arg("vector_field"),
            py::arg("properties"), py::arg("theta") = 0, py::arg("phi") = 0,
            py::arg("step_size") = 1, py::arg("do_nn_intp") = true)
       .def("set_omp_num_threads", &PliSimulator::set_omp_num_threads);
}
