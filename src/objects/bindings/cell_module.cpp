#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../cell.hpp"

namespace py = pybind11;
using Cell = object::Cell;

PYBIND11_MODULE(__cell, m) {
   py::class_<Cell>(m, "__Cell")
       // Constructors
       .def(py::init<const std::vector<float> &, const std::vector<float> &>())
       .def(py::init([](py::array_t<float, py::array::c_style> &p,
                        std::vector<float> &r) {
          auto buf_p = p.request();
          if (buf_p.ndim != 2)
             throw std::invalid_argument(
                 "array has incorrect number of dimensions: " +
                 std::to_string(buf_p.ndim) + "; expected 2");
          if (buf_p.shape[1] != 3)
             throw std::invalid_argument("array has incorrect shape: [" +
                                         std::to_string(buf_p.shape[0]) + "," +
                                         std::to_string(buf_p.shape[1]) +
                                         "]; expected [n,3]");
          if (static_cast<size_t>(buf_p.shape[0]) != r.size())
             throw std::invalid_argument("points and radii size must match");

          auto data_ptr = reinterpret_cast<float *>(buf_p.ptr);
          size_t size = std::accumulate(buf_p.shape.begin(), buf_p.shape.end(),
                                        1ULL, std::multiplies<size_t>());

          std::vector<float> points_vec(data_ptr, data_ptr + size);

          return new Cell(points_vec, r);
       }))

       // Getter
       .def_property_readonly(
           "points",
           [](const Cell &self) {
              auto data = vm::flatten(self.points());
              return py::array_t<float>(
                  std::vector<ptrdiff_t>{
                      static_cast<long long>(self.points().size()), 3},
                  data.data());
           })
       .def_property_readonly("radii",
                              [](const Cell &self) {
                                 return py::array(self.radii().size(),
                                                  self.radii().data());
                              })
       .def("size", &Cell::size)

       // Manipulators
       .def("rotate",
            [](Cell &self, const std::array<float, 9> &mat) {
               self.Rotate(vm::Mat3x3<float>(mat));
            })
       .def("rotate_around_point",
            [](Cell &self, const std::array<float, 9> &mat,
               const std::array<float, 3> &point) {
               self.RotateAroundPoint(vm::Mat3x3<float>(mat),
                                      vm::Vec3<float>(point));
            })
       .def("translate",
            [](Cell &self, const std::array<float, 3> &offset) {
               self.Translate(vm::Vec3<float>(offset));
            })
       .def("resize", &Cell::Resize)
       .def("resize_points", &Cell::ResizePoints)
       .def("resize_radii", &Cell::ResizeRadii);
}
