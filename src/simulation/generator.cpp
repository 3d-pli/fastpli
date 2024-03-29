#include "generator.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <tuple>
#include <vector>

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "fiber_bundle.hpp"
#include "include/aabb.hpp"
#include "include/omp.hpp"
#include "include/vemath.hpp"
#include "setup.hpp"

// #############################################################################
// PliGenerator Init functions
// #############################################################################

void PliGenerator::Abort(const int num) const {
   if (mpi_) {
      MPI_Abort(mpi_->comm(), num);
   } else {
      std::cerr << "ErrorCode: " << num << std::endl;
      exit(EXIT_FAILURE);
   }
}

int PliGenerator::set_omp_num_threads(int i) {
   if (i != 0) {
      if (i > omp_get_num_procs())
         omp_set_num_threads(omp_get_num_procs());
      else
         omp_set_num_threads(i);
   }

   return omp_get_max_threads();
}

void PliGenerator::SetVolume(const vm::Vec3<int64_t> global_dim,
                             const vm::Vec3<double> origin,
                             const double voxel_size) {

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
         std::cout << "rank " << mpi_->rank()
                   << ": dim.global = " << dim_.global << std::endl;
         std::cout << "rank " << mpi_->rank() << ": dim.local = " << dim_.local
                   << std::endl;
         std::cout << "rank " << mpi_->rank()
                   << ": dim.offset = " << dim_.offset << std::endl;
         std::cout << "rank " << mpi_->rank()
                   << ": dim.origin = " << dim_.origin << std::endl;
         std::cout << "rank " << mpi_->rank()
                   << ": voxel_size = " << voxel_size_ << std::endl;
      }

#ifndef NDEBUG
      if (vm::any_of(dim_.local, [&](int64_t i) { return i < 0; })) {
         std::cerr << "dim_.local" << std::endl;
         Abort(20000 + __LINE__);
      }
      if (dim_.local > dim_.global) {
         std::cerr << "dim_.global" << std::endl;
         Abort(20000 + __LINE__);
      }
#endif
   } else {
      dim_.local = global_dim;
      dim_.global = global_dim;
      dim_.origin = origin;
      dim_.offset = vm::Vec3<int64_t>(0);
      if (debug_) {
         std::cout << "dim.global = " << dim_.global << std::endl;
         std::cout << "dim.local = " << dim_.local << std::endl;
         std::cout << "dim.offset = " << dim_.offset << std::endl;
         std::cout << "dim.origin = " << dim_.origin << std::endl;
         std::cout << "voxel_size = " << voxel_size_ << std::endl;
      }
#ifndef NDEBUG
      if (dim_.local.x() < 0) {
         std::cerr << "dim_.global" << std::endl;
         Abort(20000 + __LINE__);
      }
      if (dim_.local.y() < 0) {
         std::cerr << "dim_.global" << std::endl;
         Abort(20000 + __LINE__);
      }
      if (dim_.local.z() < 0) {
         std::cerr << "dim_.global" << std::endl;
         Abort(20000 + __LINE__);
      }
      if (dim_.local > dim_.global) {
         std::cerr << "dim_.global" << std::endl;
         Abort(20000 + __LINE__);
      }
#endif
   }

   if (voxel_size <= 0)
      throw std::invalid_argument("voxel_size <= 0: " +
                                  std::to_string(voxel_size));

   voxel_size_ = voxel_size;

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

void PliGenerator::SetMPIComm(const MPI_Comm comm, const int n) {
   mpi_ = std::make_unique<MyMPI>(comm);
   if (mpi_->size() != n) {
      PyErr_WarnEx(
          PyExc_UserWarning,
          ("only " + std::to_string(mpi_->size()) + " processor found").c_str(),
          0);
   }
}

// #############################################################################
// Generation
// #############################################################################

std::tuple<std::vector<int> *, std::vector<float> *,
           std::vector<setup::PhyProps>>
PliGenerator::RunTissueGeneration(const bool only_label) {

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

   // size fibers with voxel_size
   fiber_bundles_ = fiber_bundles_org_;
   cell_populations_ = cell_populations_org_;
   if (dim_.origin != 0) {
      for (auto &fb : fiber_bundles_)
         fb.Translate(vm::cast<double>(-dim_.origin));
      for (auto &cp : cell_populations_)
         cp.Translate(vm::cast<double>(-dim_.origin));
   }

   for (auto &fb : fiber_bundles_)
      fb.Resize(1.0 / voxel_size_);
   for (auto &cp : cell_populations_)
      cp.Resize(1.0 / voxel_size_);

   if (debug_)
      std::cout << "generating fibers: " << num_fibers_ << std::endl;

#pragma omp parallel
   // every thread runs all elements, writing splits x_range
   // faster then parallel for in FillVoxelsAroundFiberSegment()
   for (size_t fb_idx = 0; fb_idx < fiber_bundles_.size(); fb_idx++) {
      const auto &fb = fiber_bundles_[fb_idx];

      if (!aabb::Overlap(volume_bb_, fb.aabb()))
         continue;

      for (size_t f_idx = 0; f_idx < fb.fibers().size(); f_idx++) {
         const auto &fiber = fb.fibers()[f_idx];

         if (fiber.size() <= 1) {
            PyErr_WarnEx(PyExc_UserWarning,
                         "fiber object with only single point", 0);
            continue;
         }

         if (!aabb::Overlap(volume_bb_, fiber.aabb()))
            continue;

         for (auto s_idx = 0u; s_idx < fiber.size() - 1; s_idx++) {
            // TODO: figure out how to incapsulate idx into fiber to only
            // parse fiber segment

            FillVoxelsAroundFiberSegment(fb_idx, f_idx, s_idx, *label_field,
                                         *vector_field, array_distance,
                                         only_label);
         }
      }
   }

   // cells
   if (debug_)
      std::cout << "generating cells: " << num_cells_ << std::endl;

#pragma omp parallel // every thread runs all elements, writing splits x_range
   for (size_t cp_idx = 0; cp_idx < cell_populations_.size(); cp_idx++) {
      const auto &cp = cell_populations_[cp_idx];

      if (!aabb::Overlap(volume_bb_, cp.aabb()))
         continue;

      for (size_t c_idx = 0; c_idx < cp.cells().size(); c_idx++) {
         const auto &cell = cp.cells()[c_idx];

         if (!aabb::Overlap(volume_bb_, cell.aabb()))
            continue;

         for (auto s_idx = 0u; s_idx < cell.size(); s_idx++) {
            // TODO: figure out how to incapsulate idx into fiber to only
            // parse fiber segment
            FillVoxelsAroundSphere(cp_idx, c_idx, s_idx, *label_field,
                                   array_distance);
         }
      }
   }

   return std::make_tuple(label_field, vector_field, GetPropertyList());
}

// #############################################################################
// Fill Cell Types
// #############################################################################

void PliGenerator::FillVoxelsAroundFiberSegment(
    const size_t fb_idx, const size_t f_idx, const size_t s_idx,
    std::vector<int> &label_field, std::vector<float> &vector_field,
    std::vector<float> &array_distance, const bool only_label) {

#ifndef NDEBUG
   if (fb_idx >= fiber_bundles_.size()) {
      std::cerr << "fb_idx" << std::endl;
      Abort(20000 + __LINE__);
   }
#endif
   const auto &fb = fiber_bundles_[fb_idx];
#ifndef NDEBUG
   if (f_idx >= fb.fibers().size()) {
      std::cerr << "f_idx" << std::endl;
      Abort(20000 + __LINE__);
   }
#endif
   const auto &fiber = fb.fibers()[f_idx];
#ifndef NDEBUG
   if (s_idx + 1 >= fiber.size()) {
      std::cerr << "s_idx" << std::endl;
      Abort(20000 + __LINE__);
   }
#endif

   const auto &p = fiber.points()[s_idx];
   const auto &q = fiber.points()[s_idx + 1];
   const auto max_radius =
       std::max(fiber.radii()[s_idx], fiber.radii()[s_idx + 1]) *
       std::sqrt(fb.layers_scale_squ().back());

   aabb::AABB<double, 3> fiber_segment_bb(p, q);
   fiber_segment_bb.min -= max_radius;
   fiber_segment_bb.max += max_radius;

   fiber_segment_bb.Intersect(volume_bb_);
   const auto min = vm::cast<int64_t>(fiber_segment_bb.min);
   const auto max = vm::cast<int64_t>(fiber_segment_bb.max + 0.5);

   double t{};
   vm::Vec3<double> min_point{};
   const auto &layers_scale_squ = fb.layers_scale_squ();

   // FIXME: TEST!!!
   for (int64_t x = std::floor(min.x()); x < max.x(); x++) {

      if (omp_in_parallel())
         // because of omp parallel for thread safe operations
         if (x % omp_get_max_threads() != omp_get_thread_num())
            continue;

      for (int64_t y = std::floor(min.y()); y < max.y(); y++) {
         for (int64_t z = std::floor(min.z()); z < max.z(); z++) {
            vm::Vec3<double> point(x + 0.5, y + 0.5,
                                   z + 0.5); // point in center of voxel -> cart

            std::tie(min_point, t) =
                ShortestPointToLineSegmentVecCalculation(point, p, q);

            const float dist_squ = vm::length2(min_point - point);
            const float r_sqr = std::pow(fiber.CalcRadius(s_idx, t), 2);
            const float &f_ly_squ = layers_scale_squ.back();
            if (dist_squ > f_ly_squ * r_sqr)
               continue;

            const size_t ind = // ind in voxel space
                (x - dim_.offset.x()) * dim_.local.y() * dim_.local.z() +
                (y - dim_.offset.y()) * dim_.local.z() + (z - dim_.offset.z());
#ifndef NDEBUG
            if (ind >= label_field.size()) {
               std::cerr << "ind >= label_field.size()" << std::endl;
               Abort(20000 + __LINE__);
            }
#endif

            if (array_distance[ind] >= dist_squ) {
               // find corresponding layer
               const auto ly_itr =
                   std::lower_bound(layers_scale_squ.begin(),
                                    layers_scale_squ.end(), dist_squ / r_sqr);
               const int ly_idx =
                   std::distance(layers_scale_squ.begin(), ly_itr);

               const int new_label = ly_idx + 1 + fb_idx * max_layer_;

               if (!flag_overlap_) {
                  if (label_field[ind] != 0) {
                     if (std::abs(label_field[ind] - new_label) >= max_layer_) {
                        flag_overlap_ = true;
                        PyErr_WarnEx(PyExc_UserWarning, "objects overlap", 0);
                     }
                  }
               }

               array_distance[ind] = dist_squ;
               label_field[ind] = new_label;

               if (!only_label) {
                  auto ly = fiber_bundles_[fb_idx].layer_orientation(ly_idx);
                  vm::Vec3<double> new_vec(0);

                  if (ly != fiber::layer::Orientation::none) {
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
#ifndef NDEBUG
   if (cp_idx >= cell_populations_.size()) {
      std::cerr << "cp_idx" << std::endl;
      Abort(20000 + __LINE__);
   }
#endif
   const auto &cp = cell_populations_[cp_idx];
#ifndef NDEBUG
   if (c_idx >= cp.cells().size()) {
      std::cerr << "c_idx" << std::endl;
      Abort(20000 + __LINE__);
   }
#endif
   const auto &cell = cp.cells()[c_idx];
#ifndef NDEBUG
   if (s_idx >= cell.size()) {
      std::cerr << "s_idx" << std::endl;
      Abort(20000 + __LINE__);
   }
#endif
   const auto &p = cell.points()[s_idx];

   aabb::AABB<double, 3> cell_sphere_bb(p - cell.radii()[s_idx],
                                        p + cell.radii()[s_idx]);
   cell_sphere_bb.Intersect(
       aabb::AABB<double, 3>(vm::cast<double>(dim_.offset),
                             vm::cast<double>(dim_.local + dim_.offset), true));
   const auto min = cell_sphere_bb.min;
   const auto max = cell_sphere_bb.max;

   const auto scale_squ = cp.scale_squ();
   for (int64_t x = std::floor(min.x()); x < std::floor(max.x()); x++) {
      // because of omp parallel
      if (omp_in_parallel())
         if (x % omp_get_max_threads() != omp_get_thread_num())
            continue;

      for (int64_t y = std::floor(min.y()); y < std::floor(max.y()); y++) {
         for (int64_t z = std::floor(min.z()); z < std::floor(max.z()); z++) {

            const vm::Vec3<double> point(x, y, z);
            const auto dist_squ = vm::length2(p - point);
            const auto r = cell.radii()[s_idx];

            if (dist_squ > scale_squ * r * r)
               continue;

            const size_t ind =
                (x - dim_.offset.x()) * dim_.local.y() * dim_.local.z() +
                (y - dim_.offset.y()) * dim_.local.z() + (z - dim_.offset.z());
#ifndef NDEBUG
            if (ind >= label_field.size()) {
               std::cerr << "ind >= label_field.size()" << std::endl;
               Abort(20000 + __LINE__);
            }
#endif

            const int new_label = cp_idx + max_layer_ * num_fiber_bundles_ +
                                  1; // +1 for background

            if (array_distance[ind] > dist_squ ||
                (array_distance[ind] == dist_squ &&
                 label_field[ind] < new_label)) {

               if (!flag_overlap_) {
                  if (label_field[ind] != 0) {
                     if (std::abs(label_field[ind] - new_label) >= max_layer_) {
                        flag_overlap_ = true;
                        PyErr_WarnEx(PyExc_UserWarning, "objects overlap", 0);
                     }
                  }
               }

               array_distance[ind] = dist_squ;
               label_field[ind] = new_label;
            }
         }
      }
   }
}

// #############################################################################
// Helper
// #############################################################################

std::tuple<vm::Vec3<double>, double>
PliGenerator::ShortestPointToLineSegmentVecCalculation(
    const vm::Vec3<double> &p, const vm::Vec3<double> &s0,
    const vm::Vec3<double> &s1) {
   const auto v = s1 - s0;
   const auto w = p - s0;

   const auto c1 = vm::dot(w, v);
   if (c1 < 0)
      return std::make_tuple(s0, 0);

   const auto c2 = vm::dot(v, v);
   if (c2 <= c1)
      return std::make_tuple(s1, 1);

   const auto b = c1 / c2;
   const auto pb = s0 + v * b;

   return std::make_tuple(pb, b);
}

std::vector<setup::PhyProps> PliGenerator::GetPropertyList() const {

   std::vector<setup::PhyProps> properties(fiber_bundles_org_.size() *
                                               max_layer_ +
                                           1 + cell_populations_org_.size());

   // background is set in __simpli.py before return simulation
   properties[0] = setup::PhyProps(0, 0); // dn has to be 0 for brackground

   for (size_t f = 0; f < fiber_bundles_org_.size(); f++) {
      for (size_t l = 0; l < fiber_bundles_org_[f].layer_size(); l++) {
         auto id = 1 + l + f * max_layer_; //+1 for brackground

         properties[id].mu = fiber_bundles_org_[f].layer_mu(l);

         if (fiber_bundles_org_[f].layer_orientation(l) ==
             fiber::layer::Orientation::none)
            properties[id].dn = 0;
         else
            properties[id].dn = fiber_bundles_org_[f].layer_dn(l);
      }
   }

   for (size_t c = 0; c < cell_populations_org_.size(); c++) {
      auto id = 1 + num_fiber_bundles_ * max_layer_ + c; //+1 for brackground
      properties[id].dn = 0;
      properties[id].mu = cell_populations_org_[c].mu();
   }

   return properties;
}
