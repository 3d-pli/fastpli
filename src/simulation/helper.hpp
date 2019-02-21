#ifndef SIMULATION_HELPER_HPP_
#define SIMULATION_HELPER_HPP_

#include "include/vemath.hpp"

struct Dimensions {
   vm::Vec3<long long> global{0};
   vm::Vec3<long long> local{0};
   vm::Vec3<long long> offset{0};
   vm::Vec3<float> origin{0};
};

struct Coordinates {
   vm::Vec3<double> tissue;
   vm::Vec2<long long> ccd;
};

#endif // SIMULATION_HELPER_HPP_
