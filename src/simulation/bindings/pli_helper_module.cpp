#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "include/vemath.hpp"
#include "pli_module.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_helper, m) {

   m.doc() = "Helper functions for simpli";

   {
      auto m_sub = m.def_submodule("rotation");
      m_sub.doc() = "Rotation Matrices";

      m_sub.def("rot_x", [](double angle) {
         auto mat = vm::rot_pi_cases::Mat3RotX<double>(angle);
         return py::array(mat.data_.size(), mat.data_.data());
      });
      m_sub.def("rot_y", [](double angle) {
         auto mat = vm::rot_pi_cases::Mat3RotY<double>(angle);
         return py::array(mat.data_.size(), mat.data_.data());
      });
      m_sub.def("rot_z", [](double angle) {
         auto mat = vm::rot_pi_cases::Mat3RotZ<double>(angle);
         return py::array(mat.data_.size(), mat.data_.data());
      });
      m_sub.def("rot_zyz", [](double a, double b, double c) {
         auto mat = vm::rot_pi_cases::Mat3RotZYZ<double>(a, b, c);
         return py::array(mat.data_.size(), mat.data_.data());
      });
   }

   {
      auto m_sub = m.def_submodule("data_container");
      m_sub.doc() = "DataContainer for simpli";

      DeclareDataContainer<int>(m_sub, "Int");
      DeclareDataContainer<float>(m_sub, "Float");
   }
}
