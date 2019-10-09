#include "generator.hpp"

#include <cassert>
#include <iostream>
#include <memory>
#include <tuple>
#include <vector>

#include "fiber_bundle.hpp"
#include "include/aabb.hpp"
#include "include/omp.hpp"
#include "include/vemath.hpp"
#include "setup.hpp"

int PliGenerator::set_omp_num_threads(int i) {
   if (i != 0) {
      if (i > omp_get_num_procs())
         omp_set_num_threads(omp_get_num_procs());
      else
         omp_set_num_threads(i);
   }

   return omp_get_max_threads();
}

void PliGenerator::SetVolume(const vm::Vec3<long long> global_dim,
                             const vm::Vec3<double> origin,
                             const double pixel_size) {

   if (global_dim.x() <= 0 || global_dim.y() <= 0 || global_dim.z() <= 0)
      throw std::invalid_argument("dim.global[any] <= 0: [" +
                                  std::to_string(global_dim.x()) + "," +
                                  std::to_string(global_dim.y()) + "," +
                                  std::to_string(global_dim.z()) + "]");

   dim_ = setup::Dimensions();
   if (mpi_) {
      mpi_->CreateCartGrid(global_dim);
      dim_ = mpi_->dim_vol();
      dim_.origin = origin;
      if (debug_) {
         std::cout << "rank " << mpi_->my_rank()
                   << ": dim.global = " << dim_.global << std::endl;
         std::cout << "rank " << mpi_->my_rank()
                   << ": dim.local = " << dim_.local << std::endl;
         std::cout << "rank " << mpi_->my_rank()
                   << ": dim.offset = " << dim_.offset << std::endl;
         std::cout << "rank " << mpi_->my_rank()
                   << ": dim.origin = " << dim_.origin << std::endl;
         std::cout << "rank " << mpi_->my_rank()
                   << ": pixel_size = " << pixel_size_ << std::endl;
      }
#ifndef NDEBUG
      if (vm::any_of(dim_.local, [&](long long i) { return i < 0; })) {
         MPI_Abort(mpi_->comm(), 2111);
      }
      if (dim_.local > dim_.global) {
         MPI_Abort(mpi_->comm(), 2112);
      }
#endif
   } else {
      dim_.local = global_dim;
      dim_.global = global_dim;
      dim_.origin = origin;
      dim_.offset = vm::Vec3<long long>(0);
      if (debug_) {
         std::cout << "dim.global = " << dim_.global << std::endl;
         std::cout << "dim.local = " << dim_.local << std::endl;
         std::cout << "dim.offset = " << dim_.offset << std::endl;
         std::cout << "dim.origin = " << dim_.origin << std::endl;
         std::cout << "pixel_size = " << pixel_size_ << std::endl;
      }
      assert(dim_.local.x() >= 0);
      assert(dim_.local.y() >= 0);
      assert(dim_.local.z() >= 0);
      assert(dim_.local <= dim_.global);
   }

   if (pixel_size <= 0)
      throw std::invalid_argument("pixel_size <= 0: " +
                                  std::to_string(pixel_size));

   pixel_size_ = pixel_size;

   volume_bb_ =
       aabb::AABB<double, 3>(vm::cast<double>(dim_.offset),
                             vm::cast<double>(dim_.local + dim_.offset), true);
}

void PliGenerator::SetFiberBundles(
    const std::vector<fiber::Bundle> &fiber_bundles) {
   fiber_bundles_org_ = fiber_bundles;

   num_fibers_ = 0;
   num_fiber_bundles_ = fiber_bundles.size();
   max_layer_ = 0;
   for (const auto &fb : fiber_bundles) {
      num_fibers_ += fb.size();
      if (static_cast<size_t>(max_layer_) < fb.layer_size())
         max_layer_ = fb.layer_size();
   }
}

void PliGenerator::SetCellPopulations(
    const std::vector<cell::Population> &cell_populations) {
   cell_populations_org_ = cell_populations;

   // num_cells_ = 0;
   num_cells_ = cell_populations.size();
   // for (const auto &cp : cell_populations)
}

void PliGenerator::SetMPIComm(const MPI_Comm comm) {
   mpi_ = std::make_unique<MyMPI>(comm);
}

std::tuple<std::vector<int> *, std::vector<float> *, setup::PhyProps>
PliGenerator::RunTissueGeneration(const bool only_label,
                                  const bool progress_bar) {

   volume_bb_ =
       aabb::AABB<double, 3>(vm::cast<double>(dim_.offset),
                             vm::cast<double>(dim_.local + dim_.offset), true);

   auto label_field = new std::vector<int>(
       dim_.local.x() * dim_.local.y() * dim_.local.z(), 0);
   auto vector_field = new std::vector<float>();

   if (!only_label)
      vector_field->resize(dim_.local.x() * dim_.local.y() * dim_.local.z() * 3,
                           0);

   // create array_distance
   std::vector<float> array_distance(dim_.local.x() * dim_.local.y() *
                                         dim_.local.z(),
                                     std::numeric_limits<float>::infinity());

   // size fibers with pixel_size
   fiber_bundles_ = fiber_bundles_org_;
   cell_populations_ = cell_populations_org_;
   if (dim_.origin != 0) {
      for (auto &fb : fiber_bundles_)
         fb.Translate(vm::cast<double>(-dim_.origin));
      for (auto &cp : cell_populations_)
         cp.Translate(vm::cast<double>(-dim_.origin));
   }

   for (auto &fb : fiber_bundles_)
      fb.Resize(1.0 / pixel_size_);
   for (auto &cp : cell_populations_)
      cp.Resize(1.0 / pixel_size_);

   int lastProgress = 0;

   if (debug_)
      std::cout << "generating fibers: " << num_fibers_ << std::endl;

   size_t progress_counter = 0;

#pragma omp parallel
   for (size_t fb_idx = 0; fb_idx < fiber_bundles_.size(); fb_idx++) {
      const auto &fb = fiber_bundles_[fb_idx];

      if (!aabb::Overlap(volume_bb_, fb.voi()))
         continue;

      for (size_t f_idx = 0; f_idx < fb.fibers().size(); f_idx++) {
         const auto &fiber = fb.fibers()[f_idx];

         if (fiber.size() <= 1)
            continue;

         if (!aabb::Overlap(volume_bb_, fiber.voi()))
            continue;

         // if (!aabb::Overlap(volume_bb_, fiber.voi()))
         //    continue;

         for (auto s_idx = 0u; s_idx < fiber.size() - 1; s_idx++) {
            // TODO: figure out how to incapsulate idx into fiber to only
            // parse fiber segment

            FillVoxelsAroundFiberSegment(fb_idx, f_idx, s_idx, *label_field,
                                         *vector_field, array_distance,
                                         only_label);
         }

         if (progress_bar) {
#pragma omp critical
            {
               progress_counter++;

               const int barWidth = 60;
               const int progress =
                   progress_counter * 100 / (num_fibers_ + num_cells_);

               if (progress - lastProgress > 1 ||
                   progress_counter == (num_fibers_ + num_cells_) ||
                   progress_counter == 1) {
                  std::cout << ": [";
                  const int pos = (barWidth * progress) / 100;
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
   }

   // cells
   if (debug_)
      std::cout << "generating cells: " << num_cells_ << std::endl;

#pragma omp parallel
   for (size_t cp_idx = 0; cp_idx < cell_populations_.size(); cp_idx++) {
      const auto &cp = cell_populations_[cp_idx];

      if (!aabb::Overlap(volume_bb_, cp.voi()))
         continue;

      for (size_t c_idx = 0; c_idx < cp.cells().size(); c_idx++) {
         const auto &cell = cp.cells()[c_idx];

         if (!aabb::Overlap(volume_bb_, cell.voi()))
            continue;

         for (auto s_idx = 0u; s_idx < cell.size(); s_idx++) {
            // TODO: figure out how to incapsulate idx into fiber to only
            // parse fiber segment
            FillVoxelsAroundSphere(cp_idx, c_idx, s_idx, *label_field,
                                   array_distance);
         }

         if (progress_bar) {
#pragma omp critical
            {
               progress_counter++;
               int barWidth = 60;
               int progress =
                   progress_counter * 100 / (num_fibers_ + num_cells_);

               if (progress - lastProgress > 1 ||
                   progress_counter == (num_fibers_ + num_cells_) ||
                   progress_counter == 1) {
                  std::cout << ": [";
                  int pos = (barWidth * progress) / 100;
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
   }

   if (progress_bar) {
      const int barWidth = 60;
      const int progress = 100;

      std::cout << ": [";
      const int pos = (barWidth * progress) / 100;
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
      std::cout << std::endl;
   }

   if (!std::any_of(label_field->begin(), label_field->end(),
                    [](int i) { return i > 0; }))
      std::cout << "WARNING: all labels are 0" << std::endl;

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

   aabb::AABB<double, 3> fiber_segment_bb(p, q);
   fiber_segment_bb.min -= max_radius;
   fiber_segment_bb.max += max_radius;

   fiber_segment_bb.Intersect(volume_bb_);
   const auto min = fiber_segment_bb.min;
   const auto max = fiber_segment_bb.max;

   double t{};
   vm::Vec3<double> min_point{};
   const auto &layers_scale_sqr = fb.layers_scale_sqr();

   for (long long x = std::floor(min.x()); x < std::floor(max.x()); x++) {

      if (omp_in_parallel())
         // because of omp parallel for thread safe operations
         if (x % omp_get_max_threads() != omp_get_thread_num())
            continue;

      for (long long y = std::floor(min.y()); y < std::floor(max.y()); y++) {
         for (long long z = std::floor(min.z()); z < std::floor(max.z()); z++) {
            vm::Vec3<double> point(x, y, z);

            std::tie(min_point, t) =
                ShortestPointToLineSegmentVecCalculation(point, p, q);

            float dist_squ = vm::length2(min_point - point);
            float ly_r = fiber.CalcRadius(s_idx, t);
            float const &f_ly_sqr = layers_scale_sqr.back();
            if (dist_squ > f_ly_sqr * ly_r * ly_r)
               continue;

            size_t ind =
                (x - dim_.offset.x()) * dim_.local.y() * dim_.local.z() +
                (y - dim_.offset.y()) * dim_.local.z() + (z - dim_.offset.z());
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
                  vm::Vec3<double> new_vec(0);

                  if (ly != fiber::layer::Orientation::background) {
                     if (ly == fiber::layer::Orientation::radial)
                        new_vec = vm::unit(point - min_point);
                     else if (ly == fiber::layer::Orientation::parallel)
                        new_vec = vm::unit(q - p);
                  }

                  std::copy(new_vec.begin(), new_vec.end(),
                            &vector_field[ind * 3]);
               }
            }
         }
      }
   }
}

void PliGenerator::FillVoxelsAroundSphere(const size_t cp_idx,
                                          const size_t c_idx,
                                          const size_t s_idx,
                                          std::vector<int> &label_field,
                                          std::vector<float> &array_distance) {

   assert(cp_idx < cell_populations_.size());
   const auto &cp = cell_populations_[cp_idx];

   assert(c_idx < cp.cells().size());
   const auto &cell = cp.cells()[c_idx];

   assert(s_idx < cell.size());
   const auto &p = cell.points()[s_idx];

   aabb::AABB<double, 3> cell_sphere_bb(p - cell.radii()[s_idx],
                                        p + cell.radii()[s_idx]);
   cell_sphere_bb.Intersect(
       aabb::AABB<double, 3>(vm::cast<double>(dim_.offset),
                             vm::cast<double>(dim_.local + dim_.offset), true));
   const auto min = cell_sphere_bb.min;
   const auto max = cell_sphere_bb.max;

   const auto scale_sqr = cp.scale_sqr();
   for (long long x = std::floor(min.x()); x < std::floor(max.x()); x++) {

      // because of omp parallel
      if (omp_in_parallel())
         if (x % omp_get_max_threads() != omp_get_thread_num())
            continue;

      for (long long y = std::floor(min.y()); y < std::floor(max.y()); y++) {
         for (long long z = std::floor(min.z()); z < std::floor(max.z()); z++) {

            vm::Vec3<double> point(x, y, z);

            auto dist_squ = vm::length2(p - point);
            auto r = cell.radii()[s_idx];

            if (dist_squ > scale_sqr * r * r)
               continue;

            size_t ind =
                (x - dim_.offset.x()) * dim_.local.y() * dim_.local.z() +
                (y - dim_.offset.y()) * dim_.local.z() + (z - dim_.offset.z());
            assert(ind < label_field.size());

            int new_label = cp_idx + max_layer_ * num_fiber_bundles_ +
                            1; // +1 for background

            if (array_distance[ind] > dist_squ ||
                (array_distance[ind] == dist_squ &&
                 label_field[ind] < new_label)) {
               array_distance[ind] = dist_squ;
               label_field[ind] = new_label;
            }
         }
      }
   }
}

std::tuple<vm::Vec3<double>, double>
PliGenerator::ShortestPointToLineSegmentVecCalculation(
    const vm::Vec3<double> &p, const vm::Vec3<double> &s0,
    const vm::Vec3<double> &s1) {
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

setup::PhyProps PliGenerator::GetPropertyList() const {

   setup::PhyProps properties(fiber_bundles_org_.size() * max_layer_ + 1 +
                              cell_populations_org_.size());

   // background
   // TODO: set background properties
   properties.dn(0) = 0;
   properties.mu(0) = 0;

   for (size_t f = 0; f < fiber_bundles_org_.size(); f++) {
      for (size_t l = 0; l < fiber_bundles_org_[f].layer_size(); l++) {
         auto id = 1 + l + f * max_layer_; //+1 for brackground

         properties.mu(id) = fiber_bundles_org_[f].layer_mu(l);

         if (fiber_bundles_org_[f].layer_orientation(l) ==
             fiber::layer::Orientation::background)
            properties.dn(id) = 0;
         else
            properties.dn(id) = fiber_bundles_org_[f].layer_dn(l);
      }
   }

   for (size_t c = 0; c < cell_populations_org_.size(); c++) {
      auto id = 1 + num_fiber_bundles_ * max_layer_ + c; //+1 for brackground
      properties.dn(id) = 0;
      properties.mu(id) = cell_populations_org_[c].mu();
   }

   return properties;
}
