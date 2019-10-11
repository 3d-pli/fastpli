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
              double voxel_size, double wavelength, double tissue_refrection,
              bool interpolate, bool untilt_sensor_view, bool flip_z,
              std::vector<double> filter_rotations) {
              setup::Setup setup;
              setup.step_size = step_size;
              setup.light_intensity = light_intensity;
              setup.voxel_size = voxel_size;
              setup.wavelength = wavelength;
              setup.tissue_refrection = tissue_refrection;
              setup.interpolate = interpolate;
              setup.untilt_sensor_view = untilt_sensor_view;
              setup.flip_z = flip_z;
              setup.filter_rotations = filter_rotations;
              self.SetSetup(setup);
           },
           py::arg("step_size"), py::arg("light_intensity"),
           py::arg("voxel_size"), py::arg("wavelength"),
           py::arg("tissue_refrection"), py::arg("interpolate"),
           py::arg("untilt_sensor_view"), py::arg("flip_z"),
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

              auto prop_vec = object::NpArray2Vector<double>(prop_array);
              std::vector<setup::PhyProps> properties;
              for (auto i = 0u; i < prop_vec.size(); i += 2)
                 properties.push_back(
                     // [2*i] -> dn, [2*i+1] -> mu
                     setup::PhyProps(prop_vec[i], prop_vec[i + 1]));

              auto image = new std::vector<double>(
                  self.RunSimulation(dim, label_container, vector_container,
                                     properties, setup::Tilting(theta, phi)));

              std::vector<size_t> dim_image =
                  vm::cast<size_t>(self.GetImageDim());

              return object::Vec2NpArray(image, dim_image);
           },
           py::arg("dim"), py::arg("label_field"), py::arg("vector_field"),
           py::arg("properties"), py::arg("theta") = 0, py::arg("phi") = 0)
       .def("set_omp_num_threads", &PliSimulator::set_omp_num_threads);
}
