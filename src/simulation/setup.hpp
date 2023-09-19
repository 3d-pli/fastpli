#ifndef SRC_SIMULATION_SETUP_HPP_
#define SRC_SIMULATION_SETUP_HPP_

#include <cassert>
#include <utility>
#include <vector>

#include "include/vemath.hpp"

namespace setup {

enum class InterpMode { nn, lerp, slerp };

struct Dimensions {
   vm::Vec3<int64_t> global{0};
   vm::Vec3<int64_t> local{0};
   vm::Vec3<int64_t> offset{0};
   vm::Vec3<double> origin{0};
};

struct Coordinates {
   vm::Vec3<double> tissue;
   vm::Vec2<int64_t> ccd;
};

struct Tilting {
   Tilting() = default;
   Tilting(double theta, double phi) : theta(theta), phi(phi) {}

   double theta{0};
   double phi{0};
};

struct PhyProps {
   PhyProps() = default;
   PhyProps(double dn, double mu) : dn(dn), mu(mu) {}

   void Check(void) const {
      if (mu < 0)
         throw std::invalid_argument("mu < 0: " + std::to_string(mu));
   }

   double dn{0};
   double mu{0};
};

struct Setup {
   void Check(void) const {
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
   }

   double step_size{1};
   double light_intensity{0};
   double voxel_size{1};
   double wavelength{0};
   double tissue_refrection{0};
   InterpMode interpolate{InterpMode::nn};
   bool untilt_sensor_view{false};
   bool flip_z{false};
   std::vector<double> filter_rotations;
};

} // namespace setup
#endif // SRC_SIMULATION_SETUP_HPP_
