#include "../fiber.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "include/aabb.hpp"
#include "include/vemath.hpp"

namespace py = pybind11;
using Fiber = data::Fiber;

PYBIND11_MODULE(_fiber_cpp, m) {
   m.doc() = "c++ fiber classes for containing fiber and fiberbundle data";

   py::class_<Fiber>(m, "_FiberCPP")
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

          return new Fiber(points_vec, r);
       }))

       // Getter
       .def_property_readonly(
           "points",
           [](const Fiber &self) {
              auto data = vm::flatten(self.points());
              return py::array_t<float>(
                  std::vector<ptrdiff_t>{
                      static_cast<long long>(self.points().size()), 3},
                  data.data());
           })
       .def_property_readonly("radii",
                              [](const Fiber &self) {
                                 return py::array(self.radii().size(),
                                                  self.radii().data());
                              })
       .def("size", &Fiber::size)

       // Manipulators
       .def("rotate",
            (void (Fiber::*)(const std::array<float, 9> &)) & Fiber::Rotate)
       .def("rotate_around_point",
            (void (Fiber::*)(const std::array<float, 9> &,
                             const std::array<float, 3> &)) &
                Fiber::RotateAroundPoint)
       .def("translate",
            (void (Fiber::*)(const std::array<float, 3> &)) & Fiber::Translate)
       .def("scale_points", &Fiber::ScalePoints)
       .def("scale_radii", &Fiber::ScaleRadii)
       .def("scale", &Fiber::Scale);
}
