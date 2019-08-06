#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <mpi.h>

#include "../generator.hpp"
#include "../helper.hpp"
#include "../setup.hpp"
#include "objects/cell.hpp"
#include "objects/fiber.hpp"
#include "objects/np_array_helper.hpp"

namespace py = pybind11;

// TODO: expose vm::vec with vec->data_ as an array

PYBIND11_MODULE(__generation, m) {
   m.doc() = "Generation of tissue for SimPLI";

   py::class_<PliGenerator>(m, "_Generator")
       .def(py::init())
       .def("set_volume",
            [](PliGenerator &self, std::array<size_t, 3> global_dim,
               std::array<float, 3> origin, double pixel_size,
               std::array<bool, 3> flip_direction) {
               self.SetVolume(vm::cast<long long>(vm::Vec3<size_t>(global_dim)),
                              origin, pixel_size, flip_direction);
            },
            py::arg("global_dim"), py::arg("origin"), py::arg("pixel_size"),
            py::arg("flip_direction") =
                std::array<bool, 3>{{false, false, false}})
       .def("set_fiber_bundles",
            [](PliGenerator &self,
               std::vector<std::vector<py::array_t<float, py::array::c_style>>>
                   fbs,
               std::vector<std::vector<fiber::layer::Property>> prs) {
               std::vector<fiber::Bundle> fiber_bundles;

               if (fbs.size() != prs.size())
                  throw py::value_error("fbs and prs not the same size");

               for (auto i = 0u; i < fbs.size(); i++) {

                  auto fiber_bundle = std::vector<std::vector<float>>();
                  for (auto f : fbs[i]) {
                     fiber_bundle.emplace_back(
                         object::NpArray2Vector<float>(f));
                  }
                  fiber_bundles.emplace_back(
                      fiber::Bundle(fiber_bundle, prs[i]));
               }

               self.SetFiberBundles(fiber_bundles);
            })
       .def("set_cell_populations",
            [](PliGenerator &self, object::CellPopulations cell_pops,
               std::vector<cell::Property> cell_prop) {
               if (cell_pops.size() != cell_prop.size())
                  throw py::value_error(
                      "cell_populations and properties not the same size");

               std::vector<cell::Population> cell_populations;
               for (size_t i = 0; i < cell_pops.size(); i++)
                  cell_populations.emplace_back(
                      cell::Population(cell_pops[i], cell_prop[i]));

               self.SetCellPopulations(cell_populations);
            })
       .def("set_mpi_comm",
            [](PliGenerator &self, long long comm_address) {
               MPI_Comm comm = *static_cast<MPI_Comm *>(
                   reinterpret_cast<void *>(comm_address));
               self.SetMPIComm(comm);
            })
       .def("run_generation",
            [](PliGenerator &self, bool only_label, bool progress_bar) {
               std::vector<int> *label_field;
               std::vector<float> *vector_field;
               setup::PhyProps properties;

               std::tie(label_field, vector_field, properties) =
                   self.RunTissueGeneration(only_label, progress_bar);

               std::vector<size_t> dim_label_field =
                   vm::cast<size_t>(self.dim_local());
               std::vector<size_t> dim_vector_field;
               if (!only_label) {
                  dim_vector_field = dim_label_field;
                  dim_vector_field.push_back(3);
               }
               std::vector<size_t> dim_prop{properties.size(), 2};

               return std::make_tuple(
                   object::Vec2NpArray(label_field, dim_label_field),
                   object::Vec2NpArray(vector_field, dim_vector_field),
                   object::Vec2NpArray(
                       new std::vector<double>(properties.vector()), dim_prop));
               std::cerr << "ZHH" << std::endl;
            },
            py::arg("only_label") = false, py::arg("progress_bar") = false)
       .def("dim_local",
            [](PliGenerator &self) { return self.dim_local().data_; })
       .def("dim_offset",
            [](PliGenerator &self) { return self.dim_offset().data_; })
       .def("set_omp_num_threads", &PliGenerator::set_omp_num_threads);

   // TODO: to objects?
   py::enum_<fiber::layer::Orientation>(m, "_Orientation")
       .value("background", fiber::layer::Orientation::background)
       .value("parallel", fiber::layer::Orientation::parallel)
       .value("radial", fiber::layer::Orientation::radial);

   py::class_<fiber::layer::Property>(m, "_LayerProperty")
       .def(py::init<float, float, float, char>())
       .def(py::init<float, float, float, fiber::layer::Orientation>())
       .def_readwrite("scale", &fiber::layer::Property::scale)
       .def_readwrite("dn", &fiber::layer::Property::dn)
       .def_readwrite("mu", &fiber::layer::Property::mu)
       .def_readwrite("layer_orientation",
                      &fiber::layer::Property::orientation);

   py::class_<cell::Property>(m, "_CellProperty").def(py::init<float, float>());
}
