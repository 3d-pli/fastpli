#ifndef SIMULATION_GENERATOR_HPP_
#define SIMULATION_GENERATOR_HPP_

#include <memory>
#include <tuple>
#include <vector>

#include "fiber_bundle.hpp"
#include "helper.hpp"
#include "include/aabb.hpp"
#include "include/vemath.hpp"
#include "my_mpi.hpp"
#include "simulator.hpp"

class PliGenerator {

 public:
   PliGenerator() = default;
   ~PliGenerator() = default;

   // setter
   void SetVolume(const vm::Vec3<long long> global_dim,
                  const vm::Vec3<float> origin, const float pixel_size,
                  const vm::Vec3<bool> flip_direction = false);
   void SetFiberBundles(const std::vector<fiber::Bundle> &fiber_bundles);

   std::tuple<std::vector<int> *, std::vector<float> *,
              std::vector<PliSimulator::PhyProp>>
   RunTissueGeneration(const bool only_label = false,
                       const bool progress_bar = false);

   // getter
   std::vector<unsigned short>
   CalcVisualLabelField(std::vector<int> label_field) const;
   Dimensions dim() { return dim_; };
   vm::Vec3<long long> dim_local() { return dim_.local; };
   vm::Vec3<long long> dim_global() { return dim_.global; };
   vm::Vec3<long long> dim_offset() { return dim_.offset; };
   vm::Vec3<float> dim_origin() { return dim_.origin; };

 private:
   const bool debug_ = true;
   double pixel_size_{0};
   Dimensions dim_;

   std::unique_ptr<MyMPI> mpi_ = std::make_unique<MyMPI>();

   vm::Vec3<bool> flip_direction_{false};

   std::vector<fiber::Bundle> fiber_bundles_org_;
   std::vector<fiber::Bundle> fiber_bundles_;

   size_t num_fibers_ = 0;
   int max_layer_ = 0;

   void FillVoxelsAroundFiberSegment(const size_t fb_idx, const size_t f_idx,
                                     const size_t s_idx,
                                     std::vector<int> &label_field,
                                     std::vector<float> &vector_field,
                                     std::vector<float> &array_distance,
                                     const bool only_label);
   std::tuple<vm::Vec3<float>, float>
   ShortestPointToLineSegmentVecCalculation(const vm::Vec3<float> &p,
                                            const vm::Vec3<float> &s0,
                                            const vm::Vec3<float> &s1);

   std::vector<PliSimulator::PhyProp> GetPropertyList() const;
};

#endif // SIMULATION_GENERATOR_HPP_
