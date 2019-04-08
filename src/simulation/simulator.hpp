#ifndef SIMULATION_SIMULATOR_HPP_
#define SIMULATION_SIMULATOR_HPP_

#include <functional>
#include <memory>
#include <vector>

#include "helper.hpp"
#include "include/omp.hpp"
#include "include/vemath.hpp"
#include "my_mpi.hpp"
#include "objects/np_array_container.hpp"

class PliSimulator {
 public:
   struct PhyProp {
      PhyProp() = default;
      PhyProp(double dn, double mu) : dn(dn), mu(mu){};

      double dn{};
      double mu{};
   };

   struct Setup {
      double light_intensity{};
      double wavelength{};
      double pixel_size{};

      bool untilt_sensor = true;
      std::vector<double> filter_rotations;

      // vm::Vec2<int> sensor_dim; // currently same dimension as x-y-volume
   };

   // PliSimulator() = default;
   PliSimulator() { set_omp_num_threads(1); }
   ~PliSimulator() = default;

   void SetPliSetup(const Setup pli_setup);
   void SetTissueProperties(const std::vector<PhyProp> &properties);

   vm::Vec3<long long> GetImageDim() {
      return vm::Vec3<long long>(
          dim_.local.x(), dim_.local.y(),
          static_cast<long long>(pli_setup_.filter_rotations.size()));
   };

   std::vector<float>
   RunSimulation(const vm::Vec3<long long> &dim,
                 object::container::NpArray<int> label_field,
                 object::container::NpArray<float> vector_field,
                 const double theta, const double phi, const double step_size,
                 const bool do_nn
                 //  , const bool flip_beam
   );

   int set_omp_num_threads(int num);

 private:
   // const bool debug_ = false;
   Setup pli_setup_{};
   Dimensions dim_{};
   object::container::NpArray<int> label_field_;
   object::container::NpArray<float> vector_field_;
   std::vector<PhyProp> properties_;

   std::vector<vm::Vec4<double>> signal_buffer_;
   // vm::Vec3<bool> flip_tissue_{false};

   std::unique_ptr<MyMPI> mpi_ = std::make_unique<MyMPI>();

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

   std::function<vm::Vec3<double>(long long, long long)>
   GetSensorToStartTransformation(const double theta, const double phi) const;

   std::vector<Coordinates> CalcPixelsUntilt(const double phi,
                                             const double theta);

   vm::Mat4x4<double> RetarderMatrix(const double beta, const double ret) const;

   bool CheckMPIHalo(const vm::Vec3<double> &local_pos,
                     const vm::Vec3<int> &shift_direct,
                     const std::vector<vm::Vec4<double>> &s_vec,
                     const Coordinates &startpos);
};

#endif // SIMULATION_SIMULATOR_HPP_
