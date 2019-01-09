#ifndef PLI_SIMULATOR_HPP_
#define PLI_SIMULATOR_HPP_

#include <functional>
#include <memory>
#include <vector>

#include "include/vemath.hpp"
#include "objects/vector_container.hpp"

struct TissueProperty {
   double dn{};
   double mu{};
};

struct PliSetup {
   double light_intensity{};
   double resolution{};
   double wavelength{};
   bool untilt_sensor = true;
   std::vector<double> filter_rotations;

   // vm::Vec2<int> sensor_dim; // currently same dimension as x-y-volume
};

class PliSimulator {

 public:
   PliSimulator() = default;
   ~PliSimulator() = default;

   void SetPliSetup(const PliSetup pli_setup);
   void SetTissue(data::VectorContainer<int> label_field_ptr,
                  data::VectorContainer<float> vector_field_ptr,
                  const std::array<int, 3> &dim,
                  const std::vector<TissueProperty> &properties,
                  const double pixel_size);
   void SetTissue(data::VectorContainer<int> label_field_ptr,
                  data::VectorContainer<float> vector_field_ptr,
                  const vm::Vec3<int> &dim,
                  const std::vector<TissueProperty> &properties,
                  const double pixel_size);

   std::vector<float> RunSimulation(const double theta, const double phi,
                                    const double step_size, const bool do_nn
                                    //  , const bool flip_beam
   );

 private:
   // bool debug_ = false;
   double pixel_size_{};
   PliSetup pli_setup_{};
   vm::Vec3<size_t> dim_{};
   data::VectorContainer<int> label_field_;
   data::VectorContainer<float> vector_field_;
   std::vector<TissueProperty> properties_;
   // vm::Vec3<bool> flip_tissue_{false};

   int GetLabel(const long long x, const long long y, const long long z) const;
   int GetLabel(const vm::Vec3<long long> p) const {
      return GetLabel(p.x(), p.y(), p.z());
   };
   int GetLabel(const double x, const double y, const double z) const {
      return GetLabel(static_cast<long long>(std::round(x)),
                      static_cast<long long>(std::round(y)),
                      static_cast<long long>(std::round(z)));
   };
   int GetLabel(const vm::Vec3<double> p) const {
      return GetLabel(p.x(), p.y(), p.z());
   };

   vm::Vec3<double> GetVec(const long long x, const long long y,
                           const long long z) const;
   vm::Vec3<double> GetVec(const vm::Vec3<long long> p) const {
      return GetVec(p.x(), p.y(), p.z());
   };
   vm::Vec3<double> GetVec(const double x, const double y, const double z,
                           const bool do_nn = true) const;
   vm::Vec3<double> GetVec(const vm::Vec3<double> p,
                           const bool do_nn = true) const {
      return GetVec(p.x(), p.y(), p.z(), do_nn);
   };

   vm::Vec3<double> InterpolateVec(const double x, const double y,
                                   const double z,
                                   const bool do_nn = true) const;
   vm::Vec3<double> InterpolateVec(const vm::Vec3<double> p,
                                   bool do_nn = true) const {
      return InterpolateVec(p.x(), p.y(), p.z(), do_nn);
   };

   vm::Vec3<double> TiltDirection(const double theta, const double phi) const;

   std::function<vm::Vec3<double>(int, int)>
   GetSensorToStartTransformation(const double theta, const double phi) const;

   vm::Mat4x4<double> RetarderMatrix(const double beta, const double ret) const;
};

#endif // PLI_SIMULATOR_HPP_
