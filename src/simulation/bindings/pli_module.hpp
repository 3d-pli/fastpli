#ifndef PLI_MODULE_HPP_
#define PLI_MODULE_HPP_

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../data_container.hpp"

// define opaque opbject to use reference
// PYBIND11_MAKE_OPAQUE(DataContainer<int>);
// PYBIND11_MAKE_OPAQUE(DataContainer<float>);
// PYBIND11_MAKE_OPAQUE(std::shared_ptr<std::vector<int>>);
// PYBIND11_MAKE_OPAQUE(std::shared_ptr<std::vector<float>>);

// TODO: consistent naming

template <typename T>
void DeclareDataContainer(pybind11::module &m, const std::string &typestr) {
   using Class = DataContainer<T>;
   std::string pyclass_name = std::string("DataContainer") + typestr;

   pybind11::class_<Class>(m, pyclass_name.c_str(), pybind11::buffer_protocol(),
                           pybind11::dynamic_attr())
       .def(pybind11::init<>())
       .def("size", [](DataContainer<T> &self) { return self.size(); })
       .def("resize",
            [](DataContainer<T> &self, size_t size) { self.resize(size); },
            pybind11::arg("size"))
       .def("fill",
            [](DataContainer<T> &self, T value) {
               std::fill(self.begin(), self.end(), value);
            },
            pybind11::arg("value"))
       .def("asarray",
            [](DataContainer<T> &self) {
               return pybind11::array(self.size(), self.data());
            })
       .def("set",
            [](DataContainer<T> &self, size_t i, T value) {
               if (i >= self.size())
                  throw pybind11::index_error("index out of bounds");
               self[i] = value;
            })
       .def("get", [](DataContainer<T> &self, size_t i) { return self[i]; })
       .def("set_array",
            [](DataContainer<T> &self, pybind11::object o) {
               self.reserve(o.ref_count());
               for (auto elm : o)
                  self.push_back(pybind11::cast<T>(elm));
            })
       .def("set_data_chunk",
            [](DataContainer<T> &self, size_t chunk_size) {
               self.SetDataChunk(chunk_size);
            })
       .def("next_data_chunk", [](DataContainer<T> &self, size_t chunk_size) {
          auto ptr = self.NextDataChunkPtr();
          if (ptr)
             return pybind11::array(chunk_size, ptr);

          return pybind11::array();
       });
}

#endif
