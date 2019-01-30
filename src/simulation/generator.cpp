#include "generator.hpp"

#include <cassert>
#include <iostream>
#include <memory>
#include <tuple>
#include <vector>

#include "fiber_class.hpp"
#include "helper.hpp"
#include "include/aabb.hpp"
#include "include/vemath.hpp"
#include "simulator.hpp"

void PliGenerator::SetVolume(const Dimensions dim, const float pixel_size,
                             const vm::Vec3<bool> flip_direction) {

   if (dim.global.x() <= 0 || dim.global.y() <= 0 || dim.global.z() <= 0)
      throw std::invalid_argument("dim.global[any] <= 0: [" +
                                  std::to_string(dim.global.x()) + "," +
                                  std::to_string(dim.global.y()) + "," +
                                  std::to_string(dim.global.z()) + "]");
   
   if (dim.local.x() <= 0 || dim.local.y() <= 0 || dim.local.z() <= 0)
      throw std::invalid_argument("dim.global[any] <= 0: [" +
                                  std::to_string(dim.local.x()) + "," +
                                  std::to_string(dim.local.y()) + "," +
                                  std::to_string(dim.local.z()) + "]");

   if (pixel_size <= 0)
      throw std::invalid_argument("pixel_size <= 0: " +
                                  std::to_string(pixel_size));

   dim_ = dim;
   pixel_size_ = pixel_size;
   flip_direction_ = flip_direction;
}
void PliGenerator::SetFiberBundles(
    const std::vector<fiber::Bundle> &fiber_bundles) {
   fiber_bundles_org_ = fiber_bundles;

   num_fibers_ = 0;
   max_layer_ = 0;
   for (const auto &fb : fiber_bundles) {
      num_fibers_ += fb.size();
      if (static_cast<size_t>(max_layer_) < fb.layer_size())
         max_layer_ = fb.layer_size();
   }
}

std::tuple<std::vector<int> *, std::vector<float> *,
           std::vector<PliSimulator::PhyProp>>
PliGenerator::RunTissueGeneration(const bool only_label,
                                  const bool progress_bar) {

   auto label_field = new std::vector<int>(dim_.x() * dim_.y() * dim_.z(), 0);
   auto vector_field = new std::vector<float>();

   if (!only_label)
      vector_field->resize(dim_.x() * dim_.y() * dim_.z() * 3, 0);

   // create array_distance
   std::vector<float> array_distance(dim_.x() * dim_.y() * dim_.z(),
                                     std::numeric_limits<float>::infinity());

   // size fibers with pixel_size
   fiber_bundles_ = fiber_bundles_org_;
   for (auto &fb : fiber_bundles_)
      fb.Resize(1.0 / pixel_size_);

   int lastProgress = 0;
   const auto volume_bb =
       aabb::AABB<float, 3>(vm::Vec3<float>(0), vm::cast<float>(dim_), true);

   size_t progress_counter = 0;

   for (size_t fb_idx = 0; fb_idx < fiber_bundles_.size(); fb_idx++) {
      const auto &fb = fiber_bundles_[fb_idx];

      for (size_t f_idx = 0; f_idx < fb.fibers().size(); f_idx++) {
         const auto &fiber = fb.fibers()[f_idx];

         progress_counter++;

         if (fiber.size() <= 1)
            continue;

         if (aabb::Overlap(volume_bb, fiber.voi())) {
            // auto volume_fiber_bb = volume_bb.Intersection(fiber.voi);

            for (auto s_idx = 0u; s_idx < fiber.size() - 1; s_idx++) {
               // TODO: figure out how to incapsulate idx into fiber to only
               // parse fiber
               FillVoxelsAroundFiberSegment(fb_idx, f_idx, s_idx, *label_field,
                                            *vector_field, array_distance,
                                            only_label);
            }
         }

         if (progress_bar) {
            int barWidth = 60;
            int progress = (progress_counter + 1) * 100 / num_fibers_;

            if (progress - lastProgress > 5 ||
                progress_counter == num_fibers_ - 1 || progress_counter == 0) {
               std::cout << ": [";
               int pos = barWidth * progress / 100;
               for (int pb = 0; pb < barWidth; ++pb) {
                  if (pb < pos)
                     std::cout << "=";
                  else if (pb == pos)
                     std::cout << ">";
                  else
                     std::cout << " ";
               }
               std::cout << "] " << progress << " %\r";
               std::cout.flush();
               lastProgress = progress;
            }
         }
      }
   }
   if (progress_bar)
      std::cout << std::endl;

   return std::make_tuple(label_field, vector_field, GetPropertyList());
}

void PliGenerator::FillVoxelsAroundFiberSegment(
    const size_t fb_idx, const size_t f_idx, const size_t s_idx,
    std::vector<int> &label_field, std::vector<float> &vector_field,
    std::vector<float> &array_distance, const bool only_label) {

   assert(fb_idx < fiber_bundles_.size());
   const auto &fb = fiber_bundles_[fb_idx];

   assert(f_idx < fb.fibers().size());
   const auto &fiber = fb.fibers()[f_idx];

   assert(s_idx + 1 < fiber.size());
   const auto &p = fiber.points()[s_idx];
   const auto &q = fiber.points()[s_idx + 1];
   const auto max_radius =
       std::max(fiber.radii()[s_idx], fiber.radii()[s_idx + 1]);

   aabb::AABB<float, 3> fiber_segment_bb(p, q);
   fiber_segment_bb.min -= max_radius;
   fiber_segment_bb.max += max_radius;
   fiber_segment_bb.Intersect(
       aabb::AABB<float, 3>(vm::Vec3<float>(0), vm::cast<float>(dim_), true));
   const auto min = fiber_segment_bb.min;
   const auto max = fiber_segment_bb.max;

   float t{};
   vm::Vec3<float> min_point{};
   const auto &layers_scale_sqr = fb.layers_scale_sqr();
   for (int x = std::round(min.x()); x < std::round(max.x()); x++) {
      for (int y = std::round(min.y()); y < std::round(max.y()); y++) {
         for (int z = std::round(min.z()); z < std::round(max.z()); z++) {
            vm::Vec3<float> point(x, y, z);

            std::tie(min_point, t) =
                ShortestPointToLineSegmentVecCalculation(point, p, q);

            auto dist_squ = vm::length2(min_point - point);
            auto ly_r = fiber.CalcRadius(s_idx, t);
            auto const &f_ly_sqr = layers_scale_sqr.back();
            if (dist_squ > f_ly_sqr * ly_r * ly_r)
               continue;

            size_t ind = x * dim_.y() * dim_.z() + y * dim_.z() + z;
            assert(ind < label_field.size());

            if (array_distance[ind] >= dist_squ) {
               // find corresponding layer
               auto ly_itr = std::lower_bound(layers_scale_sqr.begin(),
                                              layers_scale_sqr.end(),
                                              dist_squ / (ly_r * ly_r));
               int ly_idx = std::distance(layers_scale_sqr.begin(), ly_itr);

               array_distance[ind] = dist_squ;
               label_field[ind] = ly_idx + 1 + fb_idx * max_layer_;

               if (!only_label) {
                  auto ly = fiber_bundles_[fb_idx].layer_orientation(ly_idx);
                  vm::Vec3<float> tmpVec(0);

                  if (ly != fiber::layer::Orientation::background) {
                     if (ly == fiber::layer::Orientation::radial)
                        tmpVec = vm::unit(point - min_point);
                     else if (ly == fiber::layer::Orientation::parallel)
                        tmpVec = vm::unit(q - p);
                  }

                  // TODO: flip_direction for PM
                  // if (flip_direction_.x())
                  //    tmpVec.x() *= -1;

                  // if (flip_direction_.y())
                  //    tmpVec.y() *= -1;

                  // if (flip_direction_.z())
                  //    tmpVec.z() *= -1;

                  std::copy(tmpVec.begin(), tmpVec.end(),
                            &vector_field[ind * 3]);
               }
            }
         }
      }
   }
}

std::tuple<vm::Vec3<float>, float>
PliGenerator::ShortestPointToLineSegmentVecCalculation(
    const vm::Vec3<float> &p, const vm::Vec3<float> &s0,
    const vm::Vec3<float> &s1) {
   auto v = s1 - s0;
   auto w = p - s0;

   auto c1 = vm::dot(w, v);
   if (c1 < 0)
      return std::make_tuple(s0, 0);

   auto c2 = vm::dot(v, v);
   if (c2 <= c1)
      return std::make_tuple(s1, 1);

   auto b = c1 / c2;
   auto pb = s0 + v * b;

   return std::make_tuple(pb, b);
}

std::vector<unsigned short>
PliGenerator::CalcVisualLabelField(std::vector<int> label_field) const {
   // TODO: visual label_vield
   std::vector<unsigned short> vis(label_field.begin(), label_field.end());
   return vis;
}

std::vector<PliSimulator::PhyProp> PliGenerator::GetPropertyList() const {
   std::vector<PliSimulator::PhyProp> properties(
       fiber_bundles_org_.size() * max_layer_ + 1);

   // background
   // TODO: set background properties
   properties[0] = PliSimulator::PhyProp(0, 0);

   for (size_t f = 0; f < fiber_bundles_org_.size(); f++) {
      for (size_t l = 0; l < fiber_bundles_org_[f].layer_size(); l++) {
         auto dn = fiber_bundles_org_[f].layer_dn(l);
         auto mu = fiber_bundles_org_[f].layer_mu(l);

         auto id = 1 + l + f * max_layer_; //+1 for brackground
         properties.at(id) = PliSimulator::PhyProp(dn, mu);
      }
   }

   return properties;
}
