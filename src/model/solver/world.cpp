#include "world.hpp"

#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <utility>

#include <Python.h>

#include "fiber_class.hpp"
#include "include/aabb.hpp"
#include "include/omp.hpp"
#include "include/vemath.hpp"
#include "objects/fiber.hpp"
#include "oct_tree.hpp"

#if _VIS_LIBRARIES
#include "scene.hpp"
#endif //_VIS_LIBRARIES

object::FiberBundles World::get_fibers() const {
   object::FiberBundles fiber_bundles;

   if (fibers_.empty())
      return fiber_bundles;

   // push_back first element
   size_t i = 0;
   fiber_bundles.push_back(object::FiberBundle());
   fiber_bundles.back().push_back(fibers_[i]);

   // fiber order does not change, therefore map has the same order
   for (i = 1; i < fibers_.size(); i++) {
      // add fiber_bundle if required
      if (map_fb_idx_.at(i - 1).first != map_fb_idx_.at(i).first)
         fiber_bundles.push_back(object::FiberBundle());

      fiber_bundles.back().push_back(fibers_[i]);
   }

   return fiber_bundles;
}

std::vector<std::vector<std::vector<double>>> World::get_fibers_vector() const {
   std::vector<std::vector<std::vector<double>>> fiber_bundles;

   if (fibers_.empty())
      return fiber_bundles;

   // push_back first element
   size_t i = 0;
   fiber_bundles.push_back(std::vector<std::vector<double>>());
   fiber_bundles.back().push_back(fibers_[i].vector());

   // fiber order does not change, therefore map has the same order
   for (i = 1; i < fibers_.size(); i++) {
      // add fiber_bundle if required
      if (map_fb_idx_.at(i - 1).first != map_fb_idx_.at(i).first)
         fiber_bundles.push_back(std::vector<std::vector<double>>());

      fiber_bundles.back().push_back(fibers_[i].vector());
   }

   return fiber_bundles;
}

void World::set_fibers(const object::FiberBundles &fiber_bundles) {

   // free memory
   fibers_ = std::vector<geometry::Fiber>();
   map_fb_idx_ = std::map<size_t, std::pair<size_t, size_t>>();
   max_speed_ = std::numeric_limits<double>::max();

   size_t fb_idx = 0;
   for (auto const &fb : fiber_bundles) {
      size_t f_idx = 0;
      for (auto const &f : fb) {
         map_fb_idx_[fibers_.size()] = std::make_pair(fb_idx, f_idx);
         fibers_.push_back(geometry::Fiber(f, fibers_.size()));
         max_speed_ = std::min(max_speed_, fibers_.back().max_speed());
         f_idx++;
      }
      fb_idx++;
   }

   ResetObjCounter();
}

void World::set_fibers_vector(
    const std::vector<std::vector<std::vector<double>>> &fiber_bundles) {

   // free memory
   fibers_ = std::vector<geometry::Fiber>();
   map_fb_idx_ = std::map<size_t, std::pair<size_t, size_t>>();
   max_speed_ = std::numeric_limits<double>::max();

   size_t fb_idx = 0;
   for (auto const &fb : fiber_bundles) {
      size_t f_idx = 0;
      for (auto const &f : fb) {
         map_fb_idx_[fibers_.size()] = std::make_pair(fb_idx, f_idx);
         fibers_.push_back(geometry::Fiber(object::Fiber(f), fibers_.size()));
         max_speed_ = std::min(max_speed_, fibers_.back().max_speed());
         f_idx++;
      }
      fb_idx++;
   }

   ResetObjCounter();
}

void World::ResetObjCounter() {
   fiber_overlap_ = -1;
   num_col_obj_ = -1;
   num_obj_ = 0;
   max_level_ = -1;
   for (auto &f : fibers_) {
      f.set_max_speed(max_speed_);
      num_obj_ += f.ConeSize();
   }
}

int World::set_omp_num_threads(int i) {

   if (i > omp_get_num_procs())
      omp_set_num_threads(omp_get_num_procs());
   else
      omp_set_num_threads(i);

   return omp_get_max_threads();
}

bool World::ApplyBoundaryConditions(int max_steps) {
   // check fiber boundary conditions

   bool solved = fibers_.empty();
   for (; max_steps >= 0 && !solved; max_steps--) {

#pragma omp parallel for reduction(&& : solved)
      for (auto i = 0u; i < fibers_.size(); i++) {
         bool flag_length =
             fibers_[i].ApplyConeLengthConstrain(w_parameter_.obj_mean_length);
         bool flag_radius =
             fibers_[i].ApplyCurvatureConstrain(w_parameter_.obj_min_radius);
         solved = flag_length && flag_radius;
      }
   }

   num_obj_ = 0;
   for (auto const &fiber : fibers_)
      num_obj_ += fiber.ConeSize();

   return solved;
}

bool World::Step() {

   bool solved = true;
   fiber_overlap_ = 0;

   num_obj_ = 0;
   for (auto const &fiber : fibers_)
      num_obj_ += fiber.ConeSize();

#pragma omp parallel for
   // applying drag before so that velocity is an indicator for colored
   // visualization
   for (auto i = 0u; i < fibers_.size(); i++) {
      fibers_[i].Drag(w_parameter_.drag);
   }

   auto max_obj_size = 0;
   for (auto const &fiber : fibers_) {
      if (fiber.size() >= 2) {
         for (size_t i = 0; i < fiber.size() - 1; i++) {
            auto delta = vm::length(fiber.points()[i + 1] - fiber.points()[i]);
            delta += std::max(fiber.radii()[i + 1], fiber.radii()[i]);
            if (max_obj_size < delta)
               max_obj_size = delta;
         }
      }
   }

   // TODO: min_radius sollte gr;-er sein als gleichseitiges dreieck,
   // R=sqrt(3)/3.0*mean_...

   // create OctTree and calculate colliding objects in each leaf
   // TODO: num_threads
   OctTree otree(fibers_, 1.5 * max_obj_size, col_voi_);
   auto colliding_list = otree.Run();
   // std::cout << "otree.max_level(): " << otree.max_level() << std::endl;
   num_col_obj_ = colliding_list.size();

   if (!colliding_list.empty()) {
      solved = false;

      // convert set to vector for performance
      std::vector<std::array<size_t, 4>> colliding_vec(colliding_list.begin(),
                                                       colliding_list.end());

      // set speed of colliding objects
#pragma omp parallel for reduction(+ : fiber_overlap_)
      for (auto i = 0u; i < colliding_vec.size(); i++) {
         auto elm = colliding_vec[i];

         auto const [f0, f1, f2, f3, dist] =
             fibers_[elm[0]].Cone(elm[1]).PushConesApart(
                 fibers_[elm[2]].Cone(elm[3]));

         // WARNING: not thread safe
         fibers_[elm[0]].AddSpeed(elm[1], f0);
         fibers_[elm[0]].AddSpeed(elm[1] + 1, f1);
         fibers_[elm[2]].AddSpeed(elm[3], f2);
         fibers_[elm[2]].AddSpeed(elm[3] + 1, f3);

         fiber_overlap_ +=
             1 - dist / (std::max(fibers_[elm[0]].radii()[elm[1]],
                                  fibers_[elm[0]].radii()[elm[1] + 1]) +
                         std::max(fibers_[elm[2]].radii()[elm[3]],
                                  fibers_[elm[2]].radii()[elm[3] + 1]));
      }
   }

   // check fiber boundary conditions
#pragma omp parallel for reduction(&& : solved)
   for (auto i = 0u; i < fibers_.size(); i++) {
      bool flag_radius =
          fibers_[i].ApplyCurvatureConstrain(w_parameter_.obj_min_radius);
      solved = solved && flag_radius;
   }

#pragma omp parallel for
   for (auto i = 0u; i < fibers_.size(); i++)
      fibers_[i].ApplyConeLengthConstrain(w_parameter_.obj_mean_length);

   // move colliding objects
   if (!solved) {
#pragma omp parallel for
      for (auto i = 0u; i < fibers_.size(); i++) {
         fibers_[i].Move();
      }
   }

   return solved;
}

#if _VIS_LIBRARIES
void World::DrawScene() {
   if (scene_ == nullptr) {
      char arg0[] = "model.solver";
      char *argv[] = {arg0, nullptr};
      int argc = 1;
      scene_ = std::make_unique<Scene>(argc, argv);
   }
   scene_->DrawScene(fibers_);
}
#else
void World::DrawScene() {
   static bool flag = false;

   if (!flag) {
      flag = true;
      PyErr_WarnEx(PyExc_UserWarning,
                   "No OpenGl detected due build. Deactivating DrawScene()", 0);
   }
}
#endif //_VIS_LIBRARIES

/**
 * Calc Vertices for saving stl files
 */

const int kNumTubeMesh = 6;

struct Vertex {
   vm::Vec3<float> p{};
   vm::Vec3<float> n{};
};

typedef std::array<Vertex, 3> Triangle;

struct Face {
   Triangle data{};

   constexpr const Vertex &operator[](std::size_t i) const { return data[i]; }
   Vertex &operator[](std::size_t i) { return data[i]; }
};

using NSidePolygon = std::array<Vertex, kNumTubeMesh>;

size_t NumFaces(std::vector<geometry::Fiber> fibers) {
   size_t n = 0;
   for (auto const &f : fibers) {
      if (f.size() <= 1)
         continue;
      n += (f.size() - 1) * kNumTubeMesh * 2;
   }
   return n;
}

std::vector<Face> CalcTubeSkeleton(const geometry::Fiber fiber) {

   std::vector<Face> tube_faces;
   std::vector<NSidePolygon> tube_mesh;

   if (fiber.size() <= 1)
      return tube_faces;

   // calc tube mesh points
   NSidePolygon mesh_elm{};
   vm::Vec3<float> tangent_old = {0, 0, 1};
   vm::Vec3<float> p0{}, p1{}, pm{};

   for (auto i = 0u; i < fiber.size(); i++) {

      if (i == 0) {
         p0 = vm::cast<float>(fiber.points()[i]);
         p1 = vm::cast<float>(fiber.points()[i + 1]);
         pm = p0;
      } else if (i == fiber.size() - 1) {
         p0 = vm::cast<float>(fiber.points()[i - 1]);
         p1 = vm::cast<float>(fiber.points()[i]);
         pm = p1;
      } else {
         p0 = vm::cast<float>(fiber.points()[i - 1]);
         p1 = vm::cast<float>(fiber.points()[i + 1]);
         pm = vm::cast<float>(fiber.points()[i]);
      }

      auto tangent = p1 - p0;
      vm::normalize(tangent);

      // generate points
      if (i == 0) {
         // initialize circle
         for (uint k = 0; k < kNumTubeMesh; k++) {
            float t = k / (float)kNumTubeMesh;
            float theta = t * 2 * M_PI;
            mesh_elm[k].n = {cosf(theta), sinf(theta), 0.0f};
         }
      }

      // rotate old points onto new plane
      auto v = vm::cross(tangent_old, tangent);
      auto s = vm::length(v) + 1e-9;
      auto c = vm::dot(tangent_old, tangent);

      vm::Mat3x3<float> rot(
          {0, -v.z(), v.y(), v.z(), 0, -v.x(), -v.y(), v.x(), 0});

      // vm::Vec3<float>{0, -v.z, v.y};
      // vm::Vec3<float>{v.z, 0, -v.x};
      // vm::Vec3<float>{-v.y, v.x, 0};

      rot = (vm::IdentityMatrix<float, 3>() + rot) +
            (vm::dot(rot, rot) * (1 - c) / (s * s));

      for (uint k = 0; k < kNumTubeMesh; k++) {
         auto p = vm::dot(rot, mesh_elm[k].n);
         p = p / vm::length(p);

         mesh_elm[k].p = pm + p * fiber.radii()[i];
         mesh_elm[k].n = p;
      }
      tangent_old = tangent;
      tube_mesh.push_back(mesh_elm);
   }

   // calc faces
   tube_faces.reserve(fiber.size() * kNumTubeMesh * 2);
   for (auto i = 0ULL; i < tube_mesh.size() - 1; i++) {

      auto const &t0 = tube_mesh[i];
      auto const &t1 = tube_mesh[i + 1];

      Face v;
      for (auto i = 0u; i < kNumTubeMesh; i++) {

         v[0].p = t0[i % kNumTubeMesh].p;
         v[0].n = t0[i % kNumTubeMesh].n;
         v[1].p = t1[i % kNumTubeMesh].p;
         v[1].n = t1[i % kNumTubeMesh].n;
         v[2].p = t0[(i + 1) % kNumTubeMesh].p;
         v[2].n = t0[(i + 1) % kNumTubeMesh].n;
         tube_faces.push_back(v);

         v[0].p = t1[i % kNumTubeMesh].p;
         v[0].n = t1[i % kNumTubeMesh].n;
         v[1].p = t0[(i + 1) % kNumTubeMesh].p;
         v[1].n = t0[(i + 1) % kNumTubeMesh].n;
         v[2].p = t1[(i + 1) % kNumTubeMesh].p;
         v[2].n = t1[(i + 1) % kNumTubeMesh].n;
         tube_faces.push_back(v);
      }
   }

   return tube_faces;
}

void World::SaveSTL(const char *fname) {
   std::ofstream file;
   file.open(fname, std::ios::out | std::ios::binary);
   char attribute[2] = "0";

   // write header
   std::string header_info = "fastpli.model.solver";
   header_info.resize(80);
   file.write(header_info.c_str(), 80);

   // write number of faces
   auto const n = NumFaces(fibers_);
   file.write((char *)&n, 4);

   size_t n_tmp = 0;
   for (auto const &fiber : fibers_) {

      auto fiber_faces = CalcTubeSkeleton(fiber);

      for (auto const &f : fiber_faces) {

         // save normal vector
         auto normal = f[0].n + f[1].n + f[2].n;
         normal = normal / vm::length(normal);

         file.write((char *)(normal.data() + 0), 4);
         file.write((char *)(normal.data() + 1), 4);
         file.write((char *)(normal.data() + 2), 4);

         // save coordinates
         file.write((char *)(f[0].p.data() + 0), 4);
         file.write((char *)(f[0].p.data() + 1), 4);
         file.write((char *)(f[0].p.data() + 2), 4);

         file.write((char *)(f[1].p.data() + 0), 4);
         file.write((char *)(f[1].p.data() + 1), 4);
         file.write((char *)(f[1].p.data() + 2), 4);

         file.write((char *)(f[2].p.data() + 0), 4);
         file.write((char *)(f[2].p.data() + 1), 4);
         file.write((char *)(f[2].p.data() + 2), 4);

         file.write(attribute, 2);
         n_tmp++;
      }
   }

   if (n != n_tmp) {
      std::cerr << "Number of pre calc faces: " << n << std::endl;
      std::cerr << "Number of faces is wrong" << std::endl;
      exit(1);
   }

   file.close();
}
