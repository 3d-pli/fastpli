#include <pybind11/functional.h>
#include <pybind11/iostream.h>
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
              std::string interpolate, bool untilt_sensor_view, bool flip_z,
              std::vector<double> filter_rotations) {
              setup::Setup setup;
              setup.step_size = step_size;
              setup.light_intensity = light_intensity;
              setup.voxel_size = voxel_size;
              setup.wavelength = wavelength;
              setup.tissue_refrection = tissue_refrection;
              if (interpolate == "NN") {
                 setup.interpolate = setup::InterpMode::nn;
              } else if (interpolate == "Lerp") {
                 setup.interpolate = setup::InterpMode::lerp;
              } else if (interpolate == "Slerp") {
                 setup.interpolate = setup::InterpMode::slerp;
              } else {
                 throw std::invalid_argument(
                     "Only NN, Lerp or Slerp are supported");
              }
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
            [](PliSimulator &self, int64_t comm_address, int n) {
               MPI_Comm comm = *static_cast<MPI_Comm *>(
                   reinterpret_cast<void *>(comm_address));
               self.SetMPIComm(comm, n);
            })
       .def(
           "run_simulation",
           [](PliSimulator &self, std::array<int64_t, 3> dim,
              py::array_t<int, py::array::c_style> label_array,
              py::array_t<float, py::array::c_style> vector_array,
              py::array_t<double, py::array::c_style> prop_array, double theta,
              double phi) {
              py::scoped_ostream_redirect stream(
                  std::cout,                               // std::ostream&
                  py::module::import("sys").attr("stdout") // Python output
              );

              auto label_container = object::NpArray2Container(label_array);
              auto vector_container = object::NpArray2Container(vector_array);

              auto prop_vec = object::NpArray2Vector<double>(prop_array);
              std::vector<setup::PhyProps> properties;
              for (auto i = 0u; i < prop_vec.size(); i += 2)
                 properties.push_back(
                     // [2*i] -> dn, [2*i+1] -> mu
                     setup::PhyProps(prop_vec[i], prop_vec[i + 1]));

              auto image = new std::vector<double>(self.RunSimulation(
                  vm::Vec3<int64_t>(dim), label_container, vector_container,
                  properties, setup::Tilting(theta, phi)));

              std::vector<size_t> dim_image =
                  vm::cast<size_t>(self.GetImageDim());

              return object::Vec2NpArray(image, dim_image);
           },
           py::arg("dim"), py::arg("label_field"), py::arg("vector_field"),
           py::arg("properties"), py::arg("theta") = 0, py::arg("phi") = 0)
       .def("set_omp_num_threads", &PliSimulator::set_omp_num_threads)

#if _THESIS
       .def(
           "__field_interpolation",
           [](PliSimulator &self, std::array<int64_t, 3> dim,
              std::array<int64_t, 3> dim_int,
              py::array_t<int, py::array::c_style> label_array,
              py::array_t<float, py::array::c_style> vector_array,
              py::array_t<int, py::array::c_style> label_int_array,
              py::array_t<float, py::array::c_style> vector_int_array,
              std::string interpolate) {
              auto label_container = object::NpArray2Container(label_array);
              auto vector_container = object::NpArray2Container(vector_array);
              auto label_int_container =
                  object::NpArray2Container(label_int_array);
              auto vector_int_container =
                  object::NpArray2Container(vector_int_array);

              setup::InterpMode interp;
              if (interpolate == "NN") {
                 interp = setup::InterpMode::nn;
              } else if (interpolate == "Lerp") {
                 interp = setup::InterpMode::lerp;
              } else if (interpolate == "Slerp") {
                 interp = setup::InterpMode::slerp;
              } else {
                 throw std::invalid_argument(
                     "Only NN, Lerp or Slerp are supported");
              }

              py::scoped_ostream_redirect stream(
                  std::cout,                               // std::ostream&
                  py::module::import("sys").attr("stdout") // Python output
              );

              self.RunInterpolation(dim, dim_int, label_container,
                                    vector_container, label_int_container,
                                    vector_int_container, interp);
           },
           py::arg("dim"), py::arg("dim_int"), py::arg("label_field"),
           py::arg("vector_field"), py::arg("label_field_int"),
           py::arg("vector_field_int"), py::arg("interpolate"))
       .def(
           "__diff_angle",
           [](PliSimulator &self, py::array_t<float, py::array::c_style> v,
              py::array_t<float, py::array::c_style> u,
              py::array_t<float, py::array::c_style> r) {
              (void)self;
              auto v_container = object::NpArray2Container(v);
              auto u_container = object::NpArray2Container(u);
              auto r_container = object::NpArray2Container(r);

              py::scoped_ostream_redirect stream(
                  std::cout,                               // std::ostream&
                  py::module::import("sys").attr("stdout") // Python output
              );

              self.DiffAngle(v_container, u_container, r_container);
           },
           py::arg("v"), py::arg("u"), py::arg("r"))
#endif
       ;
}
