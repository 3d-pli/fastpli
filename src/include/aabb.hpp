#ifndef SRC_INCLUDE_AABB_HPP_
#define SRC_INCLUDE_AABB_HPP_

#include <algorithm>
#include <array>
#include <limits>
#include <ostream>
#include <tuple>
#include <vector>

#include "vemath.hpp"

namespace aabb {

template <typename T, int N> struct AABB {
   // constructors
   explicit AABB(const vm::Vec<T, N> &p);
   AABB(const vm::Vec<T, N> &p, const vm::Vec<T, N> &q, bool sorted = false);

   // defaults
   AABB() = default;
   AABB(AABB &&) = default;
   AABB(const AABB &) = default;
   AABB &operator=(AABB &&) = default;
   AABB &operator=(const AABB &) = default;
   ~AABB() = default;

   // data
   vm::Vec<T, N> min = vm::Vec<T, N>(std::numeric_limits<T>::infinity());
   vm::Vec<T, N> max = vm::Vec<T, N>(-std::numeric_limits<T>::infinity());

   // function
   // TODO: IsEmpty should be checked everytime?
   bool IsValid() const; // if min > max
   bool IsFinite() const;
   void Unite(const AABB<T, N> &a);
   void Unite(const vm::Vec<T, N> &v);
   void Intersect(const AABB<T, N> &a);
   void Intersect(const vm::Vec<T, N> &v);
};

// operators
template <typename T, int N>
std::ostream &operator<<(std::ostream &os, const aabb::AABB<T, N> aabb) {
   os << "[" << aabb.min << "," << aabb.max << "]";
   return os;
}

template <typename T, int N> AABB<T, N>::AABB(const vm::Vec<T, N> &p) {
   min = p;
   max = p;
}

template <typename T, int N>
AABB<T, N>::AABB(const vm::Vec<T, N> &p, const vm::Vec<T, N> &q,
                 const bool sorted) {
   if (sorted) {
      min = p;
      max = q;
   } else {
      for (int i = 0; i < N; i++) {
         if (p[i] < q[i]) {
            min[i] = p[i];
            max[i] = q[i];
         } else {
            min[i] = q[i];
            max[i] = p[i];
         }
      }
   }
}

template <typename T, int N> bool AABB<T, N>::IsValid() const {
   for (int i = 0; i < N; i++) {
      if (min[i] > max[i])
         return true;
   }
   return false;
}

template <typename T, int N> bool AABB<T, N>::IsFinite() const {
   for (int i = 0; i < N; i++) {
      if (std::isfinite(min[i]) || std::isfinite(max[i])) {
         return true;
      }
   }
   return false;
}

template <typename T, int N> void AABB<T, N>::Unite(const AABB<T, N> &a) {
   for (int i = 0; i < N; i++) {
      if (a.min[i] < min[i])
         min[i] = a.min[i];
      if (a.max[i] > max[i])
         max[i] = a.max[i];
   }
}

template <typename T, int N> void AABB<T, N>::Unite(const vm::Vec<T, N> &v) {
   for (int i = 0; i < N; i++) {
      if (v[i] < min[i])
         min[i] = v[i];
      if (v[i] > max[i])
         max[i] = v[i];
   }
}

template <typename T, int N> void AABB<T, N>::Intersect(const AABB<T, N> &a) {
   for (int i = 0; i < N; i++) {
      if (a.min[i] > min[i])
         min[i] = a.min[i];
      if (a.max[i] < max[i])
         max[i] = a.max[i];
   }
}

template <typename T, int N>
void AABB<T, N>::Intersect(const vm::Vec<T, N> &v) {
   for (int i = 0; i < N; i++) {
      if (v[i] > min[i])
         min[i] = v[i];
      if (v[i] < max[i])
         max[i] = v[i];
   }
}

// helper functions
template <typename T, int N, int M>
bool Overlap(const AABB<T, N> &a, const AABB<T, M> &b) {
   const auto n = std::min(N, M);

   for (int i = 0; i < n; i++)
      if (a.min[i] >= b.max[i] || a.max[i] <= b.min[i])
         return false;

   return true;
}

template <typename T, int N>
AABB<T, N> Union(AABB<T, N> a, const AABB<T, N> &b) {
   return a.Unite(b);
}

template <typename T, int N>
AABB<T, N> Union(AABB<T, N> a, const vm::Vec<T, N> &v) {
   return a.Unite(v);
}

template <typename T, int N>
AABB<T, N> Intersection(AABB<T, N> a, const AABB<T, N> &b) {
   return a.Intersect(b);
}

template <typename T, int N>
AABB<T, N> Intersection(AABB<T, N> a, const vm::Vec<T, N> &v) {
   return a.Intersect(v);
}

} // namespace aabb
#endif // SRC_INCLUDE_AABB_HPP_
