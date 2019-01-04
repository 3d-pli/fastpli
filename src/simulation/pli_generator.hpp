#ifndef PLI_GENERATOR_HPP_
#define PLI_GENERATOR_HPP_

#include <memory>
#include <tuple>
#include <vector>

#include "data_container.hpp"
#include "fiber_model.hpp"
#include "include/aabb.hpp"
#include "include/vemath.hpp"
#include "pli_simulator.hpp"

class PliGenerator {

 public:
   PliGenerator() = default;
   ~PliGenerator() = default;

   // setter

   void SetVolumeWithArrays(const std::array<int, 3> dim,
                            const double pixel_size,
                            const std::array<bool, 3> flip_direction =
                                std::array<bool, 3>{{false, false, false}});
   void SetVolume(const vm::Vec3<int> dim, const double pixel_size,
                  const vm::Vec3<bool> flip_direction = false);
   void SetFiberBundles(const std::vector<FiberBundle> &fiber_bundles);

   // manipulators
   void RotateFiberBundles(const std::array<double, 9> &rot_mat);
   void RotateFiberBundlesAroundPoint(const std::array<double, 9> &rot_mat,
                                      const std::array<double, 3> offset);
   void TranslateFiberBundles(const std::array<double, 3> &offset);

   std::tuple<DataContainer<int>, DataContainer<float>,
              std::vector<TissueProperty>>
   RunTissueGeneration(const bool only_label = false, const bool debug = false);

   // getter
   // std::vector<int> &GetLabelField() { return label_field_; };
   std::vector<ushort> CalcVisualLabelField(std::vector<int> label_field) const;
   // std::vector<double> &GetVectorField() { return vector_field_; };

   // double *Get() { return vector_field_.data(); };

 private:
   double pixel_size_{0};
   vm::Vec3<size_t> dim_{0};
   vm::Vec3<bool> flip_direction_{false};

   std::vector<FiberBundle> fiber_bundles_org_;
   std::vector<FiberBundle> fiber_bundles_;
   // std::vector<int> label_field_;
   // std::vector<double> vector_field_;

   size_t num_fibers_ = 0;
   int max_layer_ = 0;

   void FillVoxelsAroundFiberSegment(const size_t fb_idx, const size_t f_idx,
                                     const size_t s_idx,
                                     std::vector<int> &label_field,
                                     std::vector<float> &vector_field,
                                     std::vector<float> &array_distance,
                                     const bool only_label);
   std::tuple<vm::Vec3<double>, double>
   ShortestPointToLineSegmentVecCalculation(const vm::Vec3<double> &p,
                                            const vm::Vec3<double> &s0,
                                            const vm::Vec3<double> &s1);

   std::vector<TissueProperty> GetPropertyList() const;
};

#endif // PLI_GENERATOR_HPP_