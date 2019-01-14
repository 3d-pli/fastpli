#ifndef NP_ARRAY_HELPER_HPP_
#define NP_ARRAY_HELPER_HPP_

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <vector>

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

template <typename T> py::array Vec2NpArray(std::vector<T> *vec) {
   auto capsule = py::capsule(
       vec, [](void *vec) { delete reinterpret_cast<std::vector<T> *>(vec); });
   return py::array(vec->size(), vec->data(), capsule);
}
} // namespace object

#endif // NP_ARRAY_HELPER_HPP_
