#include "../core/fiber.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "aabb.hpp"
#include "vemath.hpp"

namespace py = pybind11;
using FiberData = object::FiberData;

PYBIND11_MODULE(_fiber_cpp, m) {
   m.doc() = "c++ fiber classes for containing fiber and fiberbundle data";

   py::class_<FiberData>(m, "_FiberCPP")
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

          return new FiberData(points_vec, r);
       }))

       // Getter
       .def_property_readonly(
           "points",
           [](const FiberData &self) {
              auto data = vm::flatten(self.points());
              return py::array_t<float>(
                  std::vector<ptrdiff_t>{
                      static_cast<long long>(self.points().size()), 3},
                  data.data());
           })
       .def_property_readonly("radii",
                              [](const FiberData &self) {
                                 return py::array(self.radii().size(),
                                                  self.radii().data());
                              })
       .def("size", &FiberData::size)
       .def("calc_radius", &FiberData::CalcRadius)

       // Manipulators
       .def("rotate_fiber",
            (void (FiberData::*)(const std::array<float, 9> &)) &
                FiberData::RotateFiber)
       .def("rotate_fiber_around_point",
            (void (FiberData::*)(const std::array<float, 9> &,
                                 const std::array<float, 3> &)) &
                FiberData::RotateFiberAroundPoint)
       .def("translate_fiber",
            (void (FiberData::*)(const std::array<float, 3> &)) &
                FiberData::TranslateFiber)
       .def("resize_fiber_pos", &FiberData::ResizeFiberPos)
       .def("resize_fiber_radii", &FiberData::ResizeFiberRadii)
       .def("resize_fiber", &FiberData::ResizeFiber);
}
