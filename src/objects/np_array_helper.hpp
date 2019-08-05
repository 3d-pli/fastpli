#ifndef NP_ARRAY_HELPER_HPP_
#define NP_ARRAY_HELPER_HPP_

#include <iostream>
#include <vector>

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "np_array_container.hpp"

/**
 * Helper function for pybind-modules to bind np_arrays and vector containers
 * into c++ and python
 */

namespace py = pybind11;
namespace object {

template <typename T>
object::container::NpArray<T>
NpArray2Container(py::array_t<T, py::array::c_style> &array) {
   auto buffer = array.request();
   assert(buffer.size >= 0);
   return object::container::NpArray<T>(
       reinterpret_cast<T *>(buffer.ptr), buffer.size,
       std::vector<size_t>(buffer.shape.begin(), buffer.shape.end()));
}

template <typename T>
std::vector<T> NpArray2Vector(py::array_t<T, py::array::c_style> &array) {

   auto buffer = array.request();
   auto size = std::accumulate(buffer.shape.begin(), buffer.shape.end(), 1ULL,
                               std::multiplies<size_t>());
   std::vector<T> data(size);

   for (auto i = 0u; i < data.size(); i++) {
      data[i] = reinterpret_cast<T *>(buffer.ptr)[i];
   }

   return data;
}

template <typename T>
py::array Vec2NpArray(std::vector<T> *vec,
                      std::vector<size_t> shape = std::vector<size_t>()) {

   if (shape.empty())
      shape.push_back(vec->size());

   assert(vec->size() == std::accumulate(shape.begin(), shape.end(), 1ULL,
                                         std::multiplies<size_t>()));

   std::vector<size_t> stride(shape.size(), 0);
   size_t elm_stride = sizeof(T);

   auto shape_it = shape.rbegin();
   auto stride_it = stride.rbegin();
   for (; stride_it != stride.rend(); stride_it++, shape_it++) {
      *stride_it = elm_stride;
      elm_stride *= *shape_it;
   }

   auto capsule = py::capsule(
       vec, [](void *vec) { delete reinterpret_cast<std::vector<T> *>(vec); });

   return py::array_t<T>(shape, stride, vec->data(), capsule);
}
} // namespace object

#endif // NP_ARRAY_HELPER_HPP_
