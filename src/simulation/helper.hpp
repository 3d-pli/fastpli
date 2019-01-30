#ifndef SIMULATION_HELPER_HPP_
#define SIMULATION_HELPER_HPP_

#include "include/vemath.hpp"

struct Dimensions {

   struct Offsets {
      vm::Vec3<long long> local{0};
      vm::Vec3<double> global{0};
   };

   vm::Vec3<long long> local{0};
   vm::Vec3<long long> global{0};
   Offsets offset;
};

#endif // SIMULATION_HELPER_HPP_
