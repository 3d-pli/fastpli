#ifndef INCLUDE_VEMATH_HPP_
#define INCLUDE_VEMATH_HPP_

#include <algorithm>
#include <array>
#define _USE_MATH_DEFINES
#include <cmath>
#include <numeric>
#include <ostream>
#include <vector>

namespace vm {

// TODO: test if int is better (e.g. negative number from user error)
template <typename T, size_t M, size_t N> struct MatMxN {
   std::array<T, M * N> data_;

   // defaults
   MatMxN() = default;
   MatMxN(MatMxN &&) = default;
   MatMxN(const MatMxN &) = default;
   MatMxN &operator=(MatMxN &&) = default;
   MatMxN &operator=(const MatMxN &) = default;
   ~MatMxN() = default;

   // constructors
   MatMxN(const T s) { data_.fill(s); }
   MatMxN(const std::array<T, M * N> &a) {
      std::copy(a.begin(), a.end(), data_.begin());
   }

   // x, y, z, w functionality
   template <typename TT>
   MatMxN(typename std::enable_if<
              std::is_convertible<TT, T>::value && M == 2 && N == 1, T>::type x,
          TT y) {
      data_[0] = x;
      data_[1] = static_cast<T>(y);
   }

   template <typename TT>
   MatMxN(typename std::enable_if<
              std::is_convertible<TT, T>::value && M == 3 && N == 1, T>::type x,
          TT y, TT z) {
      data_[0] = x;
      data_[1] = static_cast<T>(y);
      data_[2] = static_cast<T>(z);
   }

   template <typename TT>
   MatMxN(typename std::enable_if<
              std::is_convertible<TT, T>::value && M == 4 && N == 1, T>::type x,
          TT y, TT z, TT w) {
      data_[0] = x;
      data_[1] = static_cast<T>(y);
      data_[2] = static_cast<T>(z);
      data_[3] = static_cast<T>(w);
   }

   // conversion functions
   operator std::array<T, M * N>() const { return data_; }
   operator std::vector<T>() const {
      return std::vector<T>(data_.begin(), data_.end());
   }

   // iterators
   typename std::array<T, M * N>::iterator begin() { return data_.begin(); }
   typename std::array<T, M * N>::const_iterator begin() const {
      return data_.begin();
   }

   typename std::array<T, M * N>::iterator end() { return data_.end(); }
   typename std::array<T, M * N>::const_iterator end() const {
      return data_.end();
   }

   typename std::array<T, M * N>::iterator data() { return data_.data(); }
   typename std::array<T, M * N>::const_iterator data() const {
      return data_.data();
   }

   // operators
   constexpr const T &operator[](size_t i) const { return data_[i]; }
   T &operator[](size_t i) { return data_[i]; }

   template <typename U = T>
   constexpr const typename std::enable_if<N >= 2, U>::type &
   operator()(size_t i, size_t j) const {
      return data_[N * i + j];
   }

   template <typename U = T>
   typename std::enable_if<N >= 2, U>::type &operator()(size_t i, size_t j) {
      return data_[N * i + j];
   }

   MatMxN<T, M, N> operator+(const MatMxN<T, M, N> &A) const {
      MatMxN<T, M, N> result;
      std::transform(data_.begin(), data_.end(), A.begin(), result.begin(),
                     std::plus<T>());
      return result;
   }
   MatMxN<T, M, N> operator-(const MatMxN<T, M, N> &A) const {
      MatMxN<T, M, N> result;
      std::transform(data_.begin(), data_.end(), A.begin(), result.begin(),
                     std::minus<T>());
      return result;
   }
   MatMxN<T, M, N> operator*(const T f) const {
      MatMxN<T, M, N> result;
      std::transform(data_.begin(), data_.end(), result.begin(),
                     std::bind2nd(std::multiplies<T>(), f));
      return result;
   }
   MatMxN<T, M, N> operator/(const T f) const {
      MatMxN<T, M, N> result;
      std::transform(data_.begin(), data_.end(), result.begin(),
                     std::bind2nd(std::multiplies<T>(), 1 / f));
      return result;
   }
   MatMxN<T, M, N> operator-() const {
      MatMxN<T, M, N> result;
      std::transform(data_.begin(), data_.end(), result.begin(),
                     [](T elm) { return -elm; });
      return result;
   }
   MatMxN<T, M, N> &operator+=(const MatMxN<T, M, N> &u) {
      std::transform(data_.begin(), data_.end(), u.begin(), data_.begin(),
                     std::plus<T>());
      return *this;
   }
   MatMxN<T, M, N> &operator-=(const MatMxN<T, M, N> &u) {
      std::transform(data_.begin(), data_.end(), u.begin(), data_.begin(),
                     std::minus<T>());
      return *this;
   }
   MatMxN<T, M, N> &operator+=(const T f) {
      std::transform(data_.begin(), data_.end(), data_.begin(),
                     std::bind2nd(std::plus<T>(), f));
      return *this;
   }
   MatMxN<T, M, N> &operator-=(const T f) {
      std::transform(data_.begin(), data_.end(), data_.begin(),
                     std::bind2nd(std::minus<T>(), f));
      return *this;
   }
   MatMxN<T, M, N> &operator*=(const T f) {
      std::transform(data_.begin(), data_.end(), data_.begin(),
                     std::bind2nd(std::multiplies<T>(), f));
      return *this;
   }
   MatMxN<T, M, N> &operator/=(const T f) {
      std::transform(data_.begin(), data_.end(), data_.begin(),
                     std::bind2nd(std::multiplies<T>(), 1 / f));
      return *this;
   }

   bool operator==(const T a) const {
      return std::all_of(data_.begin(), data_.end(),
                         [&](T i) { return i == a; });
   }
   bool operator!=(const T a) const {
      return std::any_of(data_.begin(), data_.end(),
                         [&](T i) { return i != a; });
   }
   bool operator<=(const T a) const {
      return std::all_of(data_.begin(), data_.end(),
                         [&](T i) { return i <= a; });
   }
   bool operator<(const T a) const {
      return std::all_of(data_.begin(), data_.end(),
                         [&](T i) { return i < a; });
   }
   bool operator>=(const T a) const {
      return std::all_of(data_.begin(), data_.end(),
                         [&](T i) { return i >= a; });
   }
   bool operator>(const T a) const {
      return std::all_of(data_.begin(), data_.end(),
                         [&](T i) { return i > a; });
   }

   bool operator==(const MatMxN<T, M, N> &A) const { return data_ == A.data_; }
   bool operator!=(const MatMxN<T, M, N> &A) const { return data_ != A.data_; }
   bool operator<=(const MatMxN<T, M, N> &A) const {
      return std::equal(data_.begin(), data_.end(), A.data_.begin(),
                        [](T a, T b) -> bool { return a <= b; });
   }

   bool operator<(const MatMxN<T, M, N> &A) const {
      return std::equal(data_.begin(), data_.end(), A.data_.begin(),
                        [](T a, T b) -> bool { return a < b; });
   }

   bool operator>(const MatMxN<T, M, N> &A) const {
      return std::equal(data_.begin(), data_.end(), A.data_.begin(),
                        [](T a, T b) -> bool { return a > b; });
   }

   bool operator>=(const MatMxN<T, M, N> &A) const {
      return std::equal(data_.begin(), data_.end(), A.data_.begin(),
                        [](T a, T b) -> bool { return a >= b; });
   }

   // define accessors
   template <typename U = T>
   typename std::enable_if<M >= 1 && N == 1, U>::type &x() {
      return data_[0];
   }
   template <typename U = T>
   typename std::enable_if<M >= 1 && N == 1, U>::type const &x() const {
      return data_[0];
   }

   template <typename U = T>
   typename std::enable_if<M >= 2 && N == 1, U>::type &y() {
      return data_[1];
   }
   template <typename U = T>
   typename std::enable_if<M >= 2 && N == 1, U>::type const &y() const {
      return data_[1];
   }

   template <typename U = T>
   typename std::enable_if<M >= 3 && N == 1, U>::type &z() {
      return data_[2];
   }
   template <typename U = T>
   typename std::enable_if<M >= 3 && N == 1, U>::type const &z() const {
      return data_[2];
   }

   template <typename U = T>
   typename std::enable_if<M >= 4 && N == 1, U>::type &w() {
      return data_[3];
   }
   template <typename U = T>
   typename std::enable_if<M >= 4 && N == 1, U>::type const &w() const {
      return data_[3];
   }
};

// aliases
template <typename T, size_t N> using Vec = MatMxN<T, N, 1>;
template <typename T, size_t N> using MatNxN = MatMxN<T, N, N>;
template <typename T> using Vec2 = Vec<T, 2>;
template <typename T> using Vec3 = Vec<T, 3>;
template <typename T> using Vec4 = Vec<T, 4>;
template <typename T> using Mat2x2 = MatMxN<T, 2, 2>;
template <typename T> using Mat2x3 = MatMxN<T, 2, 3>;
template <typename T> using Mat2x4 = MatMxN<T, 2, 4>;
template <typename T> using Mat3x2 = MatMxN<T, 3, 2>;
template <typename T> using Mat3x3 = MatMxN<T, 3, 3>;
template <typename T> using Mat3x4 = MatMxN<T, 3, 4>;
template <typename T> using Mat4x2 = MatMxN<T, 4, 2>;
template <typename T> using Mat4x3 = MatMxN<T, 4, 3>;
template <typename T> using Mat4x4 = MatMxN<T, 4, 4>;

// ostream
template <typename T, size_t N>
std::ostream &operator<<(std::ostream &os, const Vec<T, N> v) {
   os << "[";
   for (size_t n = 0; n < N - 1; n++)
      os << v[n] << ",";
   os << v[N - 1] << "]";

   return os;
}

template <typename T, size_t M, size_t N>
std::ostream &operator<<(std::ostream &os, const MatMxN<T, M, N> &data) {
   if (M > 1)
      os << "[";

   for (size_t m = 0; m < M - 1; m++) {
      os << "[";
      for (size_t n = 0; n < N - 1; n++)
         os << data[N * m + n] << ",";
      os << data[N * m + (N - 1)] << "]";
   }
   os << "[";
   for (size_t n = 0; n < N - 1; n++)
      os << data[N * (M - 1) + n] << ",";
   os << data[N * (M - 1) + (N - 1)] << "]";

   if (M > 1)
      os << "]";

   return os;
}

// casting and typeing
template <typename T, typename U, size_t M, size_t N>
MatMxN<T, M, N> cast(const MatMxN<U, M, N> &u) {
   MatMxN<T, M, N> result;
   std::copy(u.begin(), u.end(), result.begin());
   return result;
}

template <typename T, size_t M, size_t N>
MatMxN<T, N, M> transpose(const MatMxN<T, M, N> &u) {
   MatMxN<T, N, M> result;
   for (size_t m = 0; m < M; m++) {
      for (size_t n = 0; n < N; n++) {
         result(n, m) = u(m, n);
      }
   }
   return result;
}

// multiplication
template <typename T, size_t N> T dot(const Vec<T, N> &v, const Vec<T, N> &u) {
   return std::inner_product(v.begin(), v.end(), u.begin(), T{});
}

template <typename T, size_t M, size_t N>
Vec<T, M> dot(const MatMxN<T, M, N> &A, const Vec<T, N> &v) {
   Vec<T, M> result;
   for (auto col = 0u; col < M; col++)
      result[col] = std::inner_product(v.begin(), v.end(), &A[col * N], T{});
   return result;
}

template <typename T, size_t M, size_t NM, size_t N>
MatMxN<T, M, N> dot(const MatMxN<T, M, NM> &A, const MatMxN<T, NM, N> &B) {
   MatMxN<T, M, N> result(0);
   for (auto m = 0u; m < M; m++)
      for (auto n = 0u; n < N; n++)
         for (auto inner = 0u; inner < NM; inner++)
            result(m, n) += A(m, inner) * B(inner, n);
   return result;
}

// elm multiplication / division
template <typename T, size_t M, size_t N>
MatMxN<T, M, N> elmmul(const MatMxN<T, M, N> &A, const MatMxN<T, M, N> &B) {
   MatMxN<T, M, N> result;
   std::transform(A.begin(), A.end(), B.begin(), result.begin(),
                  std::multiplies<T>());
   return result;
}

template <typename T, size_t M, size_t N>
MatMxN<T, M, N> elmdiv(const MatMxN<T, M, N> &A, const MatMxN<T, M, N> &B) {
   MatMxN<T, M, N> result;
   std::transform(A.begin(), A.end(), B.begin(), result.begin(),
                  std::divides<T>());
   return result;
}

// Vec functions
template <typename T> Vec<T, 3> cross(const Vec<T, 3> &u, const Vec<T, 3> &w) {
   Vec<T, 3> result;
   result[0] = u[1] * w[2] - u[2] * w[1];
   result[1] = u[2] * w[0] - u[0] * w[2];
   result[2] = u[0] * w[1] - u[1] * w[0];
   return result;
}

template <typename T, size_t N> void normalize(Vec<T, N> &u) {
   T norm2 = length2(u);
   if (norm2 > 0) {
      T invNorm = 1.0 / sqrt(norm2);
      std::for_each(u.begin(), u.end(), [&](T &elm) { elm *= invNorm; });
   }
}

template <typename T, size_t N> Vec<T, N> unit(const Vec<T, N> &u) {
   auto result = u;
   normalize(result);
   return result;
}

template <typename T, size_t N> T length2(const Vec<T, N> &u) {
   return std::inner_product(u.begin(), u.end(), u.begin(), T{});
}

template <typename T, size_t N> T length(const Vec<T, N> &u) {
   return sqrt(length2(u));
}

// additional functions
template <typename T, size_t M, size_t N> T sum(const MatMxN<T, M, N> &u) {
   return std::accumulate(u.begin(), u.end(), T(0));
}

template <typename T, size_t M, size_t N> T min(const MatMxN<T, M, N> &u) {
   return *std::min_element(u.begin(), u.end());
}

template <typename T, size_t M, size_t N> T max(const MatMxN<T, M, N> &u) {
   return *std::max_element(u.begin(), u.end());
}

template <typename T, size_t M, size_t N>
MatMxN<T, M, N> abs(MatMxN<T, M, N> u) {
   std::for_each(u.begin(), u.end(), [](T &elm) { elm = std::fabs(elm); });
   return u;
}

template <typename T, size_t M, size_t N>
MatMxN<T, M, N> round(MatMxN<T, M, N> u) {
   std::for_each(u.begin(), u.end(), [](T &elm) { elm = std::round(elm); });
   return u;
}

template <typename T, size_t M, size_t N, class UnaryFunction>
MatMxN<T, M, N> for_each(MatMxN<T, M, N> u, UnaryFunction f) {
   std::for_each(u.begin(), u.end(), f);
   return u;
}

template <typename T, size_t M, size_t N, class UnaryFunction>
bool any_of(MatMxN<T, M, N> u, UnaryFunction f) {
   return std::any_of(u.begin(), u.end(), f);
}

template <typename T, size_t M, size_t N, class UnaryFunction>
bool all_of(MatMxN<T, M, N> u, UnaryFunction f) {
   return std::all_of(u.begin(), u.end(), f);
}

template <typename T, size_t M, size_t N>
std::vector<T> flatten(const std::vector<MatMxN<T, M, N>> &data_vec) {
   std::vector<T> data(M * N * data_vec.size());
   size_t i = 0;
   for (auto const &v : data_vec) {
      std::copy(v.begin(), v.end(), data.begin() + i);
      i += M * N;
   }
   return data;
}

// special matrices
template <typename T, size_t N> MatNxN<T, N> IdentityMatrix() {
   MatNxN<T, N> U(0);
   for (auto i = 0u; i < N; i++)
      U(i, i) = 1;
   return U;
}

// Rotation operations
template <typename T> Mat3x3<T> Mat2Rot(T alpha) {
   T sin_a = sin(alpha);
   T cos_a = cos(alpha);
   return Mat2x2<T>({cos_a, -sin_a, sin_a, cos_a});
}
template <typename T> Mat3x3<T> Mat3RotX(T alpha) {
   T sin_a = sin(alpha);
   T cos_a = cos(alpha);

   return Mat3x3<T>({1, 0, 0, 0, cos_a, -sin_a, 0, sin_a, cos_a});
}
template <typename T> Mat3x3<T> Mat3RotY(T alpha) {
   T sin_a = sin(alpha);
   T cos_a = cos(alpha);
   return Mat3x3<T>({cos_a, 0, sin_a, 0, 1, 0, -sin_a, 0, cos_a});
}
template <typename T> Mat3x3<T> Mat3RotZ(T alpha) {
   T sin_a = sin(alpha);
   T cos_a = cos(alpha);
   return Mat3x3<T>({cos_a, -sin_a, 0, sin_a, cos_a, 0, 0, 0, 1});
}
template <typename T> Mat3x3<T> Mat3RotZYZ(T gamma, T beta, T alpha) {
   return dot(Mat3RotZ(gamma), dot(Mat3RotY(beta), Mat3RotZ(alpha)));
}

namespace rot_pi_cases {
template <typename T> constexpr T sin(T v) {
   if (v == 0)
      return 0;
   else if (v == M_PI_2)
      return 1;
   else if (v == M_PI)
      return 0;
   else if (v == M_PI + M_PI_2)
      return -1;
   else if (v == 2 * M_PI)
      return 0;

   return std::sin(v);
}

template <typename T> constexpr T cos(T v) {
   if (v == 0)
      return 1;
   else if (v == M_PI_2)
      return 0;
   else if (v == M_PI)
      return -1;
   else if (v == M_PI + M_PI_2)
      return 0;
   else if (v == 2 * M_PI)
      return 1;

   return std::cos(v);
}

template <typename T> Mat3x3<T> Mat2Rot(T alpha) {
   T sin_a = sin(alpha);
   T cos_a = cos(alpha);
   return Mat2x2<T>({cos_a, -sin_a, sin_a, cos_a});
}
template <typename T> Mat3x3<T> Mat3RotX(T alpha) {
   T sin_a = sin(alpha);
   T cos_a = cos(alpha);

   return Mat3x3<T>({1, 0, 0, 0, cos_a, -sin_a, 0, sin_a, cos_a});
}
template <typename T> Mat3x3<T> Mat3RotY(T alpha) {
   T sin_a = sin(alpha);
   T cos_a = cos(alpha);
   return Mat3x3<T>({cos_a, 0, sin_a, 0, 1, 0, -sin_a, 0, cos_a});
}
template <typename T> Mat3x3<T> Mat3RotZ(T alpha) {
   T sin_a = sin(alpha);
   T cos_a = cos(alpha);
   return Mat3x3<T>({cos_a, -sin_a, 0, sin_a, cos_a, 0, 0, 0, 1});
}
template <typename T> Mat3x3<T> Mat3RotZYZ(T gamma, T beta, T alpha) {
   return dot(Mat3RotZ(gamma), dot(Mat3RotY(beta), Mat3RotZ(alpha)));
}
} // namespace rot_pi_cases
} // namespace vm

#endif
