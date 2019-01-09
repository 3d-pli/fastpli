#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../vector_container.hpp"
#include "include/vemath.hpp"

namespace py = pybind11;

// TODO: try to move to contiguous np.array

template <typename T> using VectorContainer = data::VectorContainer<T>;

template <typename T>
void DeclareVectorContainer(pybind11::module &m, const std::string &typestr) {
   using Class = VectorContainer<T>;
   std::string pyclass_name = typestr;
   // std::string pyclass_name = std::string("VectorContainer") + typestr;

   pybind11::class_<Class>(m, pyclass_name.c_str(), pybind11::buffer_protocol(),
                           pybind11::dynamic_attr())
       .def(pybind11::init<>())
       .def("size", [](VectorContainer<T> &self) { return self.size(); })
       .def("resize",
            [](VectorContainer<T> &self, size_t size) { self.resize(size); },
            pybind11::arg("size"))
       .def("fill",
            [](VectorContainer<T> &self, T value) {
               std::fill(self.begin(), self.end(), value);
            },
            pybind11::arg("value"))
       .def("asarray",
            [](VectorContainer<T> &self) {
               return pybind11::array(self.size(), self.data());
            })
       .def("set",
            [](VectorContainer<T> &self, size_t i, T value) {
               if (i >= self.size())
                  throw py::index_error("index out of bounds");
               self[i] = value;
            })
       .def("get", [](VectorContainer<T> &self, size_t i) { return self[i]; })
       .def("set_array",
            [](VectorContainer<T> &self, pybind11::object o) {
               self.reserve(o.ref_count());
               for (auto elm : o)
                  self.push_back(pybind11::cast<T>(elm));
            })
       .def("set_data_chunk",
            [](VectorContainer<T> &self, size_t chunk_size) {
               self.SetDataChunk(chunk_size);
            })
       .def("next_data_chunk", [](VectorContainer<T> &self, size_t chunk_size) {
          auto ptr = self.NextDataChunkPtr();
          if (ptr)
             return pybind11::array(chunk_size, ptr);

          return pybind11::array();
       });
}

PYBIND11_MODULE(vector_container, m) {
   m.doc() = "VectorContainer for simpli";

   DeclareVectorContainer<int>(m, "Int");
   DeclareVectorContainer<float>(m, "Float");
}
