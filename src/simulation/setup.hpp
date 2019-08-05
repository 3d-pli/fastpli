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
   double light_intensity{};
   double voxel_size{};
   double wavelength{};

   std::vector<double> filter_rotations;

   bool untilt_sensor = true;
   // vm::Vec2<int> sensor_dim; // currently same dimension as x-y-volume
};

}; // namespace setup
#endif
