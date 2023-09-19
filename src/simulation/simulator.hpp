#ifndef SRC_SIMULATION_SIMULATOR_HPP_
#define SRC_SIMULATION_SIMULATOR_HPP_

#include <functional>
#include <memory>
#include <vector>

#include "include/omp.hpp"
#include "include/vemath.hpp"
#include "my_mpi.hpp"
#include "objects/np_array_container.hpp"
#include "setup.hpp"

class PliSimulator {

 public:
   // PliSimulator() = default;
   PliSimulator() { set_omp_num_threads(1); }
   ~PliSimulator() = default;

   void SetMPIComm(const MPI_Comm comm, const int n);
   void SetSetup(const setup::Setup setup);

   vm::Vec3<int64_t> GetImageDim() {
      return vm::Vec3<int64_t>(
          dim_.global.x(), dim_.global.y(),
          static_cast<int64_t>(setup_->filter_rotations.size()));
   }

   std::vector<double>
   RunSimulation(const vm::Vec3<int64_t> &dim,
                 object::container::NpArray<int> label_field,
                 object::container::NpArray<float> vector_field,
                 std::vector<setup::PhyProps> properties,
                 const setup::Tilting tilt);

   int set_omp_num_threads(int num);

#if _THESIS
   // only testing
   void RunInterpolation(const vm::Vec3<int64_t> &dim,
                         const vm::Vec3<int64_t> &dim_int,
                         object::container::NpArray<int> label_field,
                         object::container::NpArray<float> vector_field,
                         object::container::NpArray<int> label_field_int,
                         object::container::NpArray<float> vector_field_int,
                         const setup::InterpMode interpolate);
   void DiffAngle(object::container::NpArray<float> v,
                  object::container::NpArray<float> u,
                  object::container::NpArray<float> r);
#endif

 private:
#ifndef NDEBUG
   const bool debug_ = true;
#else
   const bool debug_ = false;
#endif
   const int knum_max_mpi_comm_{1000};

   setup::Dimensions dim_{};
   std::unique_ptr<const setup::Setup> setup_;
   std::unique_ptr<const std::vector<setup::PhyProps>> properties_;
   std::unique_ptr<const object::container::NpArray<int>> label_field_;
   std::unique_ptr<const object::container::NpArray<float>> vector_field_;

   std::vector<vm::Vec4<double>> stored_mpi_s_vec_;

   std::unique_ptr<MyMPI> mpi_;

   void CalculateDimensions(const vm::Vec3<int64_t> &global_dim);
   void CheckInput();

   int GetLabel(const int64_t x, const int64_t y, int64_t z) const;
   int GetLabel(const vm::Vec3<int64_t> p) const {
      return GetLabel(p.x(), p.y(), p.z());
   }
   int GetLabel(const double x, const double y, const double z) const;
   int GetLabel(const vm::Vec3<double> p) const {
      return GetLabel(p.x(), p.y(), p.z());
   }

   vm::Vec3<double> GetVec(const int64_t x, const int64_t y, int64_t z) const;
   vm::Vec3<double> GetVec(const vm::Vec3<int64_t> p) const {
      return GetVec(p.x(), p.y(), p.z());
   }
   vm::Vec3<double> GetVec(const double x, const double y, const double z,
                           const setup::InterpMode interpolate) const;
   vm::Vec3<double> GetVec(const vm::Vec3<double> p,
                           const setup::InterpMode interpolate) const {
      return GetVec(p.x(), p.y(), p.z(), interpolate);
   }
   vm::Vec3<double> LerpVec(const double x, const double y,
                            const double z) const;
   vm::Vec3<double> SlerpVec(const double x, const double y,
                             const double z) const;

   vm::Vec3<double> LightDirectionUnitVector(const setup::Tilting tilt) const;
   vm::Vec3<int>
   LightDirectionComponent(const vm::Vec3<double> &direction_vec) const;
   std::vector<setup::Coordinates>
   CalcStartingLightPositions(const setup::Tilting &tilt);
   std::vector<setup::Coordinates>
   CalcStartingLightPositionsTilted(const setup::Tilting &tilt);
   std::vector<setup::Coordinates>
   CalcStartingLightPositionsUntilted(const setup::Tilting &tilt);

   bool CheckMPIHalo(const vm::Vec3<double> &local_pos,
                     const vm::Vec2<int64_t> &ccd_pos,
                     const vm::Vec3<int> &shift_direct,
                     const std::vector<vm::Vec4<double>> &s_vec);

   void Abort(const int num) const;
};

#endif // SRC_SIMULATION_SIMULATOR_HPP_
