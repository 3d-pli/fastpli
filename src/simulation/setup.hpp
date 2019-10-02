#ifndef SIMULATION_SETUP_HPP_
#define SIMULATION_SETUP_HPP_

#include <cassert>
#include <utility>
#include <vector>

#include "include/vemath.hpp"

namespace setup {

struct Dimensions {
   vm::Vec3<long long> global{0};
   vm::Vec3<long long> local{0};
   vm::Vec3<long long> offset{0};
   vm::Vec3<double> origin{0};
};

struct Coordinates {
   vm::Vec3<double> tissue;
   vm::Vec2<long long> ccd;
};

struct Tilting {
   Tilting(const double theta, const double phi) : theta(theta), phi(phi) {
      if (std::abs(theta) >= M_PI_2)
         throw std::invalid_argument("abs(theta) >= pi/2: " +
                                     std::to_string(theta));
   };

   // FIXME: Tilting should be const, not theta and phi
   const double theta;
   const double phi;
};

class PhyProps {
 public:
   PhyProps() = default;
   PhyProps(size_t i) { data_.resize(2 * i); };
   PhyProps(std::vector<double> data) : data_(data){};

   std::pair<double, double> operator[](size_t i) {
      assert(2 * i + 1 < data_.size());
      return std::make_pair(dn(i), mu(i));
   }

   double &dn(size_t i) {
      assert(2 * i + 0 < data_.size());
      return data_[2 * i + 0];
   }
   double &mu(size_t i) {
      assert(2 * i + 1 < data_.size());
      return data_[2 * i + 1];
   }
   size_t size() { return data_.size() / 2; }
   void push_back(double dn, double mu) {
      data_.push_back(dn);
      data_.push_back(mu);
   };

   std::vector<double> vector() { return data_; }

 private:
   std::vector<double> data_;
};

struct Setup {

   Setup(const double step_size, const double light_intensity,
         const double voxel_size, const double wavelength,
         const double tissue_refrection, const bool interpolate,
         const bool untilt_sensor_view, const bool flip_z,
         const std::vector<double> filter_rotations)
       : step_size(step_size), light_intensity(light_intensity),
         voxel_size(voxel_size), wavelength(wavelength),
         tissue_refrection(tissue_refrection), interpolate(interpolate),
         untilt_sensor_view(untilt_sensor_view), flip_z(flip_z),
         filter_rotations(filter_rotations) {

      if (step_size <= 0)
         throw std::invalid_argument("step_size <= 0: " +
                                     std::to_string(step_size));

      if (light_intensity < 0)
         throw std::invalid_argument("light intensity < 0: " +
                                     std::to_string(light_intensity));

      if (voxel_size <= 0)
         throw std::invalid_argument("voxel_size <= 0: " +
                                     std::to_string(voxel_size));

      if (wavelength <= 0)
         throw std::invalid_argument("wavelength <= 0: " +
                                     std::to_string(wavelength));

      if (tissue_refrection < 1)
         throw std::invalid_argument("tissue_refrection < 1: " +
                                     std::to_string(tissue_refrection));

      if (filter_rotations.empty())
         throw std::invalid_argument("filter_rotations is empty: []");
   };

   const double step_size;
   const double light_intensity;
   const double voxel_size;
   const double wavelength;
   const double tissue_refrection;
   const bool interpolate;
   const bool untilt_sensor_view;
   const bool flip_z;
   const std::vector<double> filter_rotations;
};

} // namespace setup
#endif
