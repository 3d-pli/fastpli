#ifndef VECTOR_CONTAINER_HPP_
#define VECTOR_CONTAINER_HPP_

#include <exception>
#include <memory>
#include <vector>

namespace object {
namespace container {
template <typename T> struct Vector {
   Vector() = default;
   ~Vector() = default;
   Vector(std::vector<T> &v) { (*data_) = v; }
   Vector(std::shared_ptr<std::vector<T>> &sv) { data_ = sv; }

   // data
   std::shared_ptr<std::vector<T>> data_ = std::make_shared<std::vector<T>>();
   size_t chunk_size_ = 0;
   size_t i_chunk_ = 0;

   // operators
   constexpr const T &operator[](std::size_t i) const {
      return data_->operator[](i);
   }
   T &operator[](std::size_t i) { return data_->operator[](i); }

   // getter
   std::vector<T> ptr() { return *data_; }

   // vector functions
   void clear() { data_->clear(); }
   T *data() { return data_->data(); }
   void resize(size_t i) { data_->resize(i); }
   void reserve(size_t i) { data_->reserve(i); }
   void push_back(T val) { return data_->push_back(val); }
   size_t size() const { return data_->size(); }

   // iterators
   typename std::vector<T>::iterator begin() { return data_->begin(); }
   typename std::vector<T>::const_iterator begin() const {
      return data_.begin();
   }

   typename std::vector<T>::iterator end() { return data_->end(); }
   typename std::vector<T>::const_iterator end() const { return data_.end(); }

   // chunking data for h5io
   void SetDataChunk(size_t chunk_size);
   T *NextDataChunkPtr();
};

template <typename T> void Vector<T>::SetDataChunk(size_t chunk_size) {
   if (data_->size() % chunk_size != 0)
      throw std::invalid_argument("chunk size must be dividable by data size");

   chunk_size_ = chunk_size;
   i_chunk_ = 0;
}

template <typename T> T *Vector<T>::NextDataChunkPtr() {
   auto ptr = data_->data() + chunk_size_ * i_chunk_;
   if (chunk_size_ * i_chunk_ >= data_->size())
      ptr = nullptr;

   i_chunk_++;
   return ptr;
}
} // namespace container
} // namespace object

#endif
