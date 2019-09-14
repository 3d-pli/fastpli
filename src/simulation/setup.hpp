#ifndef SIMULATION_SETUP_HPP_
#define SIMULATION_SETUP_HPP_

#include <cassert>
#include <utility>
#include <vector>

namespace setup {
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
         const bool interpolation, const bool untilt_sensor,
         const std::vector<double> filter_rotations)
       : step_size(step_size), light_intensity(light_intensity),
         voxel_size(voxel_size), wavelength(wavelength),
         interpolation(interpolation), untilt_sensor(untilt_sensor),
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

      if (filter_rotations.empty())
         throw std::invalid_argument("filter_rotations is empty: []");
   };

   const double step_size;
   const double light_intensity;
   const double voxel_size;
   const double wavelength;
   const bool interpolation;
   const bool untilt_sensor;
   const std::vector<double> filter_rotations;
};

} // namespace setup
#endif
