#ifndef SIMULATION_GENERATOR_HPP_
#define SIMULATION_GENERATOR_HPP_

#include <memory>
#include <tuple>
#include <vector>

#include "cell_population.hpp"
#include "fiber_bundle.hpp"
#include "include/aabb.hpp"
#include "include/omp.hpp"
#include "include/vemath.hpp"
#include "my_mpi.hpp"
#include "setup.hpp"

class PliGenerator {

 public:
   // PliGenerator() = default;
   PliGenerator() { set_omp_num_threads(1); }
   ~PliGenerator() = default;

   // setter
   void SetVolume(const vm::Vec3<long long> global_dim,
                  const vm::Vec3<double> origin, const double pixel_size);
   void SetFiberBundles(const std::vector<fiber::Bundle> &fiber_bundles);
   void
   SetCellPopulations(const std::vector<cell::Population> &cell_populations);
   void SetMPIComm(const MPI_Comm comm);

   std::tuple<std::vector<int> *, std::vector<float> *, setup::PhyProps>
   RunTissueGeneration(const bool only_label = false,
                       const bool progress_bar = false);
   int set_omp_num_threads(int num);

   // getter
   setup::Dimensions dim() const { return dim_; };
   vm::Vec3<long long> dim_local() const { return dim_.local; };
   vm::Vec3<long long> dim_global() const { return dim_.global; };
   vm::Vec3<long long> dim_offset() const { return dim_.offset; };
   vm::Vec3<double> dim_origin() const { return dim_.origin; };

 private:
#ifndef NDEBUG
   const bool debug_ = true;
#else
   const bool debug_ = false;
#endif

   double pixel_size_{0};
   setup::Dimensions dim_;
   aabb::AABB<double, 3> volume_bb_{};

   std::unique_ptr<MyMPI> mpi_ = nullptr;

   std::vector<fiber::Bundle> fiber_bundles_org_;
   std::vector<fiber::Bundle> fiber_bundles_;

   std::vector<cell::Population> cell_populations_org_;
   std::vector<cell::Population> cell_populations_;

   size_t num_fibers_ = 0;
   size_t num_fiber_bundles_ = 0;
   int max_layer_ = 0;
   size_t num_cells_ = 0;

   void FillVoxelsAroundFiberSegment(const size_t fb_idx, const size_t f_idx,
                                     const size_t s_idx,
                                     std::vector<int> &label_field,
                                     std::vector<float> &vector_field,
                                     std::vector<float> &array_distance,
                                     const bool only_label);
   void FillVoxelsAroundSphere(const size_t cp_idx, const size_t c_idx,
                               const size_t s_idx,
                               std::vector<int> &label_field,
                               std::vector<float> &array_distance);
   std::tuple<vm::Vec3<double>, double>
   ShortestPointToLineSegmentVecCalculation(const vm::Vec3<double> &p,
                                            const vm::Vec3<double> &s0,
                                            const vm::Vec3<double> &s1);

   setup::PhyProps GetPropertyList() const;
};

#endif // SIMULATION_GENERATOR_HPP_
