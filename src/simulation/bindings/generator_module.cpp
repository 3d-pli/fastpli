#include <exception>
#include <iostream>
#include <tuple>
#include <vector>

#include <pybind11/functional.h>
#include <pybind11/iostream.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <mpi.h>

#include "../generator.hpp"
#include "../setup.hpp"
#include "objects/cell.hpp"
#include "objects/fiber.hpp"
#include "objects/np_array_helper.hpp"

namespace py = pybind11;

PYBIND11_MODULE(__generation, m) {
   m.doc() = "Generation of tissue for SimPLI";

   py::class_<PliGenerator>(m, "_Generator")
       .def(py::init())
       .def(
           "set_volume",
           [](PliGenerator &self, std::array<size_t, 3> global_dim,
              std::array<double, 3> origin, double voxel_size) {
              self.SetVolume(vm::cast<long long>(vm::Vec3<size_t>(global_dim)),
                             origin, voxel_size);
           },
           py::arg("global_dim"), py::arg("origin"), py::arg("voxel_size"))
       .def(
           "set_fiber_bundles",
           [](PliGenerator &self,
              std::vector<std::vector<py::array_t<double, py::array::c_style>>>
                  fbs,
              std::vector<std::vector<std::tuple<double, double, double, char>>>
                  prs) {
              std::vector<fiber::Bundle> fiber_bundles;

              if (fbs.size() != prs.size())
                 throw py::value_error("fbs and prs not the same size");

              for (auto i = 0u; i < fbs.size(); i++) {
                 auto fiber_bundle = std::vector<std::vector<double>>();
                 for (auto f : fbs[i]) {
                    fiber_bundle.emplace_back(
                        object::NpArray2Vector<double>(f));
                 }

                 std::vector<fiber::layer::Property> prs_elm;
                 for (auto const &p : prs[i])
                    prs_elm.push_back(
                        fiber::layer::Property(std::get<0>(p), std::get<1>(p),
                                               std::get<2>(p), std::get<3>(p)));

                 fiber_bundles.emplace_back(
                     fiber::Bundle(fiber_bundle, prs_elm));
              }

              self.SetFiberBundles(fiber_bundles);
           })
       .def("set_cell_populations",
            [](PliGenerator &self,
               std::vector<std::vector<py::array_t<double, py::array::c_style>>>
                   cell_pops,
               std::vector<std::tuple<double, double>> cell_prop) {
               std::vector<cell::Population> cell_populations;
               if (cell_pops.size() != cell_prop.size())
                  throw py::value_error(
                      "cell_populations and properties not the same size");

               for (auto i = 0u; i < cell_pops.size(); i++) {
                  auto cell_population = std::vector<std::vector<double>>();
                  for (auto c : cell_pops[i]) {
                     cell_population.emplace_back(
                         object::NpArray2Vector<double>(c));
                  }

                  cell_populations.emplace_back(cell::Population(
                      cell_population,
                      cell::Property(std::get<0>(cell_prop[i]),
                                     std::get<1>(cell_prop[i]))));
               }

               self.SetCellPopulations(cell_populations);
            })
       .def("set_mpi_comm",
            [](PliGenerator &self, long long comm_address, int n) {
               MPI_Comm comm = *static_cast<MPI_Comm *>(
                   reinterpret_cast<void *>(comm_address));
               self.SetMPIComm(comm, n);
            })
       .def(
           "run_generation",
           [](PliGenerator &self, bool only_label) {
              py::scoped_ostream_redirect stream(
                  std::cout,                               // std::ostream&
                  py::module::import("sys").attr("stdout") // Python output
              );

              std::vector<int> *label_field;
              std::vector<float> *vector_field;
              std::vector<setup::PhyProps> properties;

              std::tie(label_field, vector_field, properties) =
                  self.RunTissueGeneration(only_label);

              std::vector<size_t> dim_label_field =
                  vm::cast<size_t>(self.dim_local());
              std::vector<size_t> dim_vector_field;
              if (!only_label) {
                 dim_vector_field = dim_label_field;
                 dim_vector_field.push_back(3);
              }
              std::vector<size_t> dim_prop{properties.size(), 2};

              std::vector<double> prop_vec;
              for (auto elm : properties) {
                 prop_vec.push_back(elm.dn);
                 prop_vec.push_back(elm.mu);
              }

              return std::make_tuple(
                  object::Vec2NpArray(label_field, dim_label_field),
                  object::Vec2NpArray(vector_field, dim_vector_field),
                  object::Vec2NpArray(new std::vector<double>(prop_vec),
                                      dim_prop));
           },
           py::arg("only_label") = false)
       .def("dim_local",
            [](PliGenerator &self) { return self.dim_local().data_; })
       .def("dim_offset",
            [](PliGenerator &self) { return self.dim_offset().data_; })
       .def("set_omp_num_threads", &PliGenerator::set_omp_num_threads);
}
