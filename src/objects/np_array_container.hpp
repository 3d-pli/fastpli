#ifndef SRC_OBJECTS_NP_ARRAY_CONTAINER_HPP_
#define SRC_OBJECTS_NP_ARRAY_CONTAINER_HPP_

// #include <exception>
#include <cassert>
#include <functional>
#include <iostream>
#include <vector>

namespace object {
namespace container {
template <typename T> class NpArray {
   /**
    * This Container ist for reading a numpy array. Therefore python will
    * dealocate the memory. Be sure, that the data is passed at runtime of the
    * function, since python can change the memory adress anytime.
    */
 public:
   NpArray(T *ptr, size_t n, std::vector<size_t> shape)
       : data_{ptr}, size_(n), shape_(shape) {

      assert(shape.size() > 0);
      assert(std::accumulate(shape.begin(), shape.end(), 1ULL,
                             std::multiplies<size_t>()) == n);
   }

   NpArray() = default;

   // operators
   constexpr const T &operator[](size_t i) const {
      assert(i < size_);
      return data_[i];
   }
   T &operator[](size_t i) {
      assert(i < size_);
      return data_[i];
   }

   // getter
   size_t size() const { return size_; }
   size_t ndim() const { return shape_.size(); }
   const std::vector<size_t> &shape() const { return shape_; }

 private:
   T *const data_;
   size_t size_ = 0;
   std::vector<size_t> shape_{};
};

} // namespace container
} // namespace object

#endif // SRC_OBJECTS_NP_ARRAY_CONTAINER_HPP_
