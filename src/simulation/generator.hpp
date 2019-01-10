#ifndef SIMULATION_GENERATOR_HPP_
#define SIMULATION_GENERATOR_HPP_

#include <memory>
#include <tuple>
#include <vector>

#include "fiber_class.hpp"
#include "include/aabb.hpp"
#include "include/vemath.hpp"
#include "objects/vector_container.hpp"
#include "simulator.hpp"

class PliGenerator {

 public:
   PliGenerator() = default;
   ~PliGenerator() = default;

   // setter

   void SetVolumeWithArrays(const std::array<int, 3> dim,
                            const float pixel_size,
                            const std::array<bool, 3> flip_direction =
                                std::array<bool, 3>{{false, false, false}});
   void SetVolume(const vm::Vec3<int> dim, const float pixel_size,
                  const vm::Vec3<bool> flip_direction = false);
   void SetFiberBundles(const std::vector<fiber::Bundle> &fiber_bundles);

   std::tuple<object::container::Vector<int>, object::container::Vector<float>,
              std::vector<PliSimulator::TissueProperty>>
   RunTissueGeneration(const bool only_label = false,
                       const bool progress_bar = false);

   // getter
   std::vector<unsigned short>
   CalcVisualLabelField(std::vector<int> label_field) const;

 private:
   double pixel_size_{0};
   vm::Vec3<size_t> dim_{0};
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

   std::vector<PliSimulator::TissueProperty> GetPropertyList() const;
};

#endif // SIMULATION_GENERATOR_HPP_
