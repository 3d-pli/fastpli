#include "simulator.hpp"

#include <cmath>
#include <exception>
#include <functional>
#include <iostream>
#include <memory>
#include <numeric>
#include <vector>

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "include/omp.hpp"
#include "include/vemath.hpp"
#include "my_mpi.hpp"
#include "objects/np_array_container.hpp"
#include "setup.hpp"

namespace pi_case = vm::rot_pi_cases;

// #############################################################################
// Optical Elements
// #############################################################################

vm::Mat4x4<double> PolX(const double p) {
   // see dissertation hendrik wiese
   vm::Mat4x4<double> M = {{1, p, 0, 0, p, 1, 0, 0, 0, 0, sqrt(1 - p * p), 0, 0,
                            0, 0, sqrt(1 - p * p)}};
   return M * 0.5;
}

vm::Mat4x4<double> PolY(const double p) {
   // see dissertation hendrik wiese
   vm::Mat4x4<double> M = {{1, -p, 0, 0, -p, 1, 0, 0, 0, 0, sqrt(1 - p * p), 0,
                            0, 0, 0, sqrt(1 - p * p)}};
   return M * 0.5;
}

vm::Mat4x4<double> RetarderMatrix(const double beta, const double ret) {
   vm::Mat4x4<double> M = {
       {1, 0, 0, 0, 0,
        pow(pi_case::cos(beta), 2) +
            pi_case::cos(ret) * pow(pi_case::sin(beta), 2),
        (1 - pi_case::cos(ret)) * pi_case::sin(beta) * pi_case::cos(beta),
        -pi_case::sin(ret) * pi_case::sin(beta), 0,
        (1 - pi_case::cos(ret)) * pi_case::sin(beta) * pi_case::cos(beta),
        pow(pi_case::sin(beta), 2) +
            pi_case::cos(ret) * pow(pi_case::cos(beta), 2),
        pi_case::sin(ret) * pi_case::cos(beta), 0,
        pi_case::sin(ret) * pi_case::sin(beta),
        -pi_case::sin(ret) * pi_case::cos(beta), pi_case::cos(ret)}};
   return M;
}

// #############################################################################
// PliSimulator Init functions
// #############################################################################

void PliSimulator::Abort(const int num) const {
   if (mpi_) {
      MPI_Abort(mpi_->comm(), num);
   } else {
      std::cerr << "ErrorCode: " << num << std::endl;
      exit(EXIT_FAILURE);
   }
}

int PliSimulator::set_omp_num_threads(int i) {
   if (i != 0) {
      if (i > omp_get_num_procs())
         omp_set_num_threads(omp_get_num_procs());
      else
         omp_set_num_threads(i);
   }
   return omp_get_max_threads();
}

void PliSimulator::SetMPIComm(const MPI_Comm comm) {
   auto mpi = std::make_unique<MyMPI>(comm);
   if (mpi->size() > 1)
      mpi_ = std::move(mpi);
   else
      PyErr_WarnEx(PyExc_UserWarning, "only one processor found", 0);
}

void PliSimulator::SetSetup(const setup::Setup setup) {
   setup_ = std::make_unique<setup::Setup>(setup);
}

void PliSimulator::CalculateDimensions(const vm::Vec3<long long> &global_dim) {

   if (global_dim.x() <= 0 || global_dim.y() <= 0 || global_dim.z() <= 0) {
      std::cerr << "global_dim[any] <= 0: [" + std::to_string(global_dim.x()) +
                       "," + std::to_string(global_dim.y()) + "," +
                       std::to_string(global_dim.z()) + "]"
                << std::endl;
      Abort(30000 + __LINE__);
   }

   if (mpi_) {
      mpi_->CreateCartGrid(global_dim);
      dim_ = mpi_->dim_vol();
      mpi_->set_num_rho(setup_->filter_rotations.size());
   } else {
      dim_.local = global_dim;
      dim_.global = global_dim;
      dim_.offset = vm::Vec3<long long>(0);
   }

   if (dim_.local.x() * dim_.local.y() * dim_.local.z() * 1ULL !=
       label_field_->size()) {
      std::cerr << "dim_.local.x() * dim_.local.y() * "
                   "dim_.local.z() != label_field_.size()"
                << std::endl;
      Abort(30000 + __LINE__);
   }

   if (dim_.local.x() * dim_.local.y() * dim_.local.z() * 3ULL !=
       vector_field_->size()) {
      std::cerr << "dim_.local.x() * dim_.local.y() * "
                   "dim_.local.z()* 3 != vector_field_.size()"
                << std::endl;
      Abort(30000 + __LINE__);
   }
}

void PliSimulator::CheckInput() {

   if (!setup_) {
      std::cerr << "pli_setup not set yet" << std::endl;
      Abort(30000 + __LINE__);
   }
   setup_->Check();

   if (label_field_->size() == 0) {
      std::cerr << "label_field.size() == 0" << std::endl;
      Abort(30000 + __LINE__);
   }

   if (label_field_->size() * 3 != vector_field_->size()) {
      std::cerr << "label_field_.size()*3 != vector_field_.size()" << std::endl;
      Abort(30000 + __LINE__);
   }

   int min = std::numeric_limits<int>::max();
   int max = std::numeric_limits<int>::min();

// warning: can take a significant amount of time
#pragma omp parallel for reduction(min : min) reduction(max : max)
   for (size_t i = 0; i < label_field_->size(); i++) {
      min = std::min(min, (*label_field_)[i]);
      max = std::max(max, (*label_field_)[i]);
   }

   if (min < 0 || max < 0) {
      std::cerr << "label < 0 detected" << std::endl;
      Abort(30000 + __LINE__);
   }

   if (static_cast<size_t>(max) >= properties_->size()) {
      std::cerr << "label exceed properties.size()" << std::endl;
      Abort(30000 + __LINE__);
   }

   if ((*properties_)[0].dn != 0) {
      std::cerr << "background birefringence has to be 0" << std::endl;
      Abort(30000 + __LINE__);
   }

   for (auto const &elm : *properties_)
      elm.Check();
}

// #############################################################################
// Simulation
// #############################################################################

std::vector<double>
PliSimulator::RunSimulation(const vm::Vec3<long long> &global_dim,
                            object::container::NpArray<int> label_field,
                            object::container::NpArray<float> vector_field,
                            std::vector<setup::PhyProps> properties,
                            const setup::Tilting tilt) {

   // transfer data to class
   label_field_ =
       std::make_unique<object::container::NpArray<int>>(label_field);
   vector_field_ =
       std::make_unique<object::container::NpArray<float>>(vector_field);
   properties_ = std::make_unique<std::vector<setup::PhyProps>>(properties);

   CheckInput();
   CalculateDimensions(global_dim);

   const auto n_rho = setup_->filter_rotations.size();
   const double lambda = setup_->wavelength * 1e-9;
   const double thickness = setup_->voxel_size * 1e-6 * setup_->step_size;

   // polarizer and lambda/4 retarder
   const double polarization_x = 1; // TODO: via setup_
   const double polarization_y = 1; // TODO: via setup_

   const auto polarizer_x = PolX(polarization_x);
   const auto polarizer_y = PolY(polarization_y);
   const auto m_lambda_4 = RetarderMatrix(M_PI_2, -M_PI_2);

   // initial values
   vm::Vec4<double> signal_0 = {{setup_->light_intensity, 0, 0, 0}};
   signal_0 = vm::dot(m_lambda_4, vm::dot(polarizer_x, signal_0));
   const std::vector<vm::Vec4<double>> s_vec_0(n_rho, signal_0);

   const vm::Mat3x3<double> rotation =
       vm::rot_pi_cases::Mat3RotZYZ(tilt.phi, tilt.theta, -tilt.phi);

   const vm::Vec3<double> light_dir_vec = LightDirectionUnitVector(tilt);
   const vm::Vec3<int> light_dir_comp = LightDirectionComponent(light_dir_vec);
   const vm::Vec3<double> light_step = light_dir_vec * setup_->step_size;

   std::vector<double> intensity_signal(
       dim_.global.x() * dim_.global.y() * n_rho,
       -std::numeric_limits<double>::infinity());

   auto light_positions = CalcStartingLightPositions(tilt);
   // add half step so that the coordinate is in the center of its current light
   // path
   for (auto &lp : light_positions)
      lp.tissue += light_step * 0.5;

   int num_total_mpi_comm = 0;

   while (!light_positions.empty()) {
#pragma omp parallel for schedule(guided)
      for (size_t s = 0; s < light_positions.size(); s++) {

         bool flag_save_ray = true;

         auto local_pos =
             light_positions[s].tissue - vm::cast<double>(dim_.offset);
         const auto ccd_pos = light_positions[s].ccd;
         auto s_vec = s_vec_0;

         if (!stored_mpi_s_vec_.empty()) {
#ifndef NDEBUG
            if ((s + 1) * n_rho > stored_mpi_s_vec_.size())
               Abort(30000 + __LINE__);
#endif
            std::copy(stored_mpi_s_vec_.begin() + s * n_rho,
                      stored_mpi_s_vec_.begin() + (s + 1) * n_rho,
                      s_vec.begin());
         }

         // go inside loop as long midpoint of current step is inside tissue.
         // local_pos.x and local_pos.y are guaranteed to be safe by
         // CalcStartingLightPositions()
         for (; std::floor(local_pos.z()) >= 0 &&
                std::floor(local_pos.z()) < dim_.local.z();
              local_pos += light_step) {

            // check if communication is neccesary
            if (CheckMPIHalo(local_pos, ccd_pos, light_dir_comp, s_vec)) {
               flag_save_ray = false;
               break;
            }

            const auto label = GetLabel(local_pos);
#ifndef NDEBUG
            if (label >= static_cast<long long>(properties_->size()))
               Abort(30000 + __LINE__);
#endif

            // calculate physical parameters
            const auto dn = (*properties_)[label].dn;
            const auto mu = (*properties_)[label].mu * 1e3;
            const auto attenuation = pow(exp(-0.5 * mu * thickness), 2);

            if (dn == 0 || label == 0) {
               // label == 0 if outside tissue -> no vector
               if (mu == 0)
                  continue;

               for (auto rho = 0u; rho < n_rho; rho++)
                  s_vec[rho] *= attenuation;
               continue;
            }

            auto vec = GetVec(local_pos, setup_->interpolate);

            if (vec == 0) {
               PyErr_WarnEx(PyExc_UserWarning,
                            "optical_axis is 0 for tissue > 0", 0);
               continue;
            }
#ifndef NDEBUG
            if (std::isnan(vec.x()) || std::isnan(vec.y()) ||
                std::isnan(vec.z()))
               Abort(30000 + __LINE__);
#endif

            const auto d_rel = dn * 4.0 * thickness / lambda * vm::length(vec);

            if (tilt.theta != 0)
               vec = vm::dot(rotation, vec);

            // calculate spherical coordinates
            const auto alpha = asin(vec.z() / vm::length(vec));
            const auto ret = M_PI_2 * d_rel * pow(cos(alpha), 2.0);
            const auto sin_r = sin(ret);
            const auto cos_r = cos(ret);

            const auto phii = atan2(vec.y(), vec.x());

            for (auto rho = 0u; rho < n_rho; rho++) {
               const auto beta = 2 * (setup_->filter_rotations[rho] - phii);
               const auto sin_b = sin(-beta);
               const auto cos_b = cos(-beta);

               const auto a1 = s_vec[rho][1] * cos_b - s_vec[rho][2] * sin_b;
               const auto c1 = s_vec[rho][1] * sin_b + s_vec[rho][2] * cos_b;
               const auto b1 = c1 * cos_r + s_vec[rho][3] * sin_r;

               // is equivalent to (R*M*R*S)*att
               s_vec[rho] = {{s_vec[rho][0], a1 * cos_b + b1 * sin_b,
                              -a1 * sin_b + b1 * cos_b,
                              -c1 * sin_r + s_vec[rho][3] * cos_r}};
               s_vec[rho] *= attenuation;

#ifndef NDEBUG
               if (std::isnan(s_vec[rho][0]) || std::isnan(s_vec[rho][1]) ||
                   std::isnan(s_vec[rho][2]) || std::isnan(s_vec[rho][3]))
                  Abort(30000 + __LINE__);
#endif
            }
         }

         if (flag_save_ray) {
            // save only, if light ray has reached the end of the volume
            size_t ccd_idx = ccd_pos.x() * dim_.global.y() + ccd_pos.y();

#ifndef NDEBUG
            if (ccd_idx * n_rho >= intensity_signal.size())
               Abort(30000 + __LINE__);
#endif
            for (auto rho = 0u; rho < n_rho; rho++) {
               s_vec[rho] = vm::dot(polarizer_y, s_vec[rho]);
#ifndef NDEBUG
               if (std::isnan(s_vec[rho][0]))
                  Abort(30000 + __LINE__);
#endif
               intensity_signal[ccd_idx * n_rho + rho] = s_vec[rho][0];
            }
         }
      }

      light_positions.clear();

      if (mpi_) {
         // stay in communication loop until all processes do not have to
         // communicate anymore.
         int num_remaining_light_elm = -1;
         while (num_remaining_light_elm != 0 && light_positions.empty()) {
            std::tie(light_positions, stored_mpi_s_vec_) =
                mpi_->CommunicateData();
            num_remaining_light_elm =
                mpi_->AllreduceSum(light_positions.size());
            if (debug_) {
               std::cout << "rank " << mpi_->rank()
                         << ": num_remaining_light_elm: "
                         << num_remaining_light_elm << std::endl;
            }
            num_total_mpi_comm++;
            if (num_total_mpi_comm > knum_max_mpi_comm_) {
               std::cerr << "mpi: num_total_mpi_comm > "
                            "num_max_total_communications."
                         << ": light_positions.size(): "
                         << light_positions.size() << std::endl;
               Abort(30000 + __LINE__);
            }
         }
      }
   }
#ifndef NDEBUG
   if (mpi_) {
      std::vector<int> mask(intensity_signal.size() / n_rho);
      for (size_t i = 0; i < mask.size(); i++)
         mask[i] = std::isnormal(intensity_signal[i * n_rho]);
      mask = mpi_->AllreduceSum(mask);

      for (size_t i = 0; i < mask.size(); i++) {
         if (mask[i] != 1) {
            Abort(30000 + __LINE__);
         }
      }
   }
#endif

   if (mpi_)
      intensity_signal = mpi_->AllreduceMax(intensity_signal);

   return intensity_signal;
}

// #############################################################################
// Get Label
// #############################################################################

auto llfloor = [](double x) { return static_cast<long long>(std::floor(x)); };
auto llceil = [](double x) { return static_cast<long long>(std::ceil(x)); };

int PliSimulator::GetLabel(const double x, const double y,
                           const double z) const {
   return GetLabel(llfloor(x), llfloor(y), llfloor(z));
}

int PliSimulator::GetLabel(const long long x, const long long y,
                           long long z) const {

   if (setup_->flip_z)
      z = dim_.local.z() - 1 - z;

   if (x < 0 || x >= dim_.local.x() || y < 0 || y >= dim_.local.y() || z < 0 ||
       z >= dim_.local.z())
      return 0; // Outside of tissue -> background

   const size_t idx =
       x * dim_.local.y() * dim_.local.z() + y * dim_.local.z() + z;

#ifndef NDEBUG
   if (idx >= label_field_->size())
      Abort(30000 + __LINE__);
#endif

   return (*label_field_)[idx];
}

// #############################################################################
// Get Vector
// #############################################################################

vm::Vec3<double> VectorOrientationAddition(const vm::Vec3<double> v,
                                           const vm::Vec3<double> u) {

   if (vm::dot(v, u) >= 0)
      return v + u;
   else
      return v - u;
}

vm::Vec3<double> LinearOrientationInterpolation(const vm::Vec3<double> v,
                                                const vm::Vec3<double> u,
                                                const double t) {
   assert(t >= 0 && t <= 1);
   assert(t >= 0 && t <= 1);
   assert(t >= 0 && t <= 1);

   if (v == 0)
      return u;
   if (u == 0)
      return v;

   return VectorOrientationAddition(v * (1 - t), u * t);
}

vm::Vec3<double> SphericalInterpolation(vm::Vec3<double> v, vm::Vec3<double> u,
                                        const double t) {

   assert(t >= 0 && t <= 1);

   if (v == 0)
      return u;
   if (u == 0)
      return v;
   // vm::normalize(v);
   // vm::normalize(u);
   if (vm::dot(v, u) < 0)
      u = -u;
   if (vm::length2(v - u) < 1e-12)
      return v;

   const auto theta = acos(vm::dot(v, u) / (vm::length(v) * vm::length(u)));
   assert(std::isfinite(theta));

   return v * sin((1 - t) * theta) / sin(theta) +
          u * sin(t * theta) / sin(theta);
}

vm::Vec3<double>
PliSimulator::GetVec(const double x, const double y, const double z,
                     const setup::InterpMode interpolate) const {
   // Vec[i,j,k] is at position [i+0.5, j+0.5, k+0.5]

   if (interpolate == setup::InterpMode::nn)
      return GetVec(llfloor(x), llfloor(y), llfloor(z));

   const auto label =
       GetLabel(llfloor(x), llfloor(y), llfloor(z)); // nearest neighbor
   if (label == 0) {
      return vm::Vec3<double>(0, 0, 0);
   }

   const auto x0 =
       std::min(std::max(llfloor(x - 0.5), 0LL), dim_.local.x() - 1);
   const auto y0 =
       std::min(std::max(llfloor(y - 0.5), 0LL), dim_.local.y() - 1);
   const auto z0 =
       std::min(std::max(llfloor(z - 0.5), 0LL), dim_.local.z() - 1);

   const auto x1 = std::min(std::max(llceil(x - 0.5), 0LL), dim_.local.x() - 1);
   const auto y1 = std::min(std::max(llceil(y - 0.5), 0LL), dim_.local.y() - 1);
   const auto z1 = std::min(std::max(llceil(z - 0.5), 0LL), dim_.local.z() - 1);

   if (x0 == x1 && y0 == y1 && z0 == z1)
      return GetVec(x0, y0, z0);

   auto xd = (x - x0 - 0.5) / (x1 - x0); // TODO: dont need x1-x0
   auto yd = (y - y0 - 0.5) / (y1 - y0);
   auto zd = (z - z0 - 0.5) / (z1 - z0);

   if (x0 == x1)
      xd = 0;
   if (y0 == y1)
      yd = 0;
   if (z0 == z1)
      zd = 0;

   // only interpolate if same label id
   const auto c000 = GetVec(x0, y0, z0) * (GetLabel(x0, y0, z0) == label);
   const auto c100 = GetVec(x1, y0, z0) * (GetLabel(x1, y0, z0) == label);
   const auto c010 = GetVec(x0, y1, z0) * (GetLabel(x0, y1, z0) == label);
   const auto c110 = GetVec(x1, y1, z0) * (GetLabel(x1, y1, z0) == label);
   const auto c001 = GetVec(x0, y0, z1) * (GetLabel(x0, y0, z1) == label);
   const auto c101 = GetVec(x1, y0, z1) * (GetLabel(x1, y0, z1) == label);
   const auto c011 = GetVec(x0, y1, z1) * (GetLabel(x0, y1, z1) == label);
   const auto c111 = GetVec(x1, y1, z1) * (GetLabel(x1, y1, z1) == label);

   auto interp = &LinearOrientationInterpolation;
   if (interpolate == setup::InterpMode::slerp)
      interp = &SphericalInterpolation;

   const auto c00 = interp(c000, c100, xd);
   const auto c01 = interp(c001, c101, xd);
   const auto c10 = interp(c010, c110, xd);
   const auto c11 = interp(c011, c111, xd);

   const auto c0 = interp(c00, c10, yd);
   const auto c1 = interp(c01, c11, yd);

   return interp(c0, c1, zd);
}

vm::Vec3<double> PliSimulator::GetVec(const long long x, const long long y,
                                      long long z) const {
   if (setup_->flip_z)
      z = dim_.local.z() - 1 - z;

   const size_t idx =
       x * dim_.local.y() * dim_.local.z() * 3 + y * dim_.local.z() * 3 + z * 3;

#ifndef NDEBUG
   assert(x >= 0);
   assert(y >= 0);
   assert(z >= 0);
   assert(x < dim_.local.x());
   assert(y < dim_.local.y());
   assert(z < dim_.local.z());
   if (idx >= vector_field_->size()) {
      std::cerr << "idx >= vector_field_->size()" << std::endl;
      Abort(30000 + __LINE__);
   }
#endif
   return vm::Vec3<double>((*vector_field_)[idx], (*vector_field_)[idx + 1],
                           (*vector_field_)[idx + 2]);
}

// #############################################################################
// Light Path functions
// #############################################################################

vm::Vec3<double>
PliSimulator::LightDirectionUnitVector(const setup::Tilting tilt) const {
   return vm::Vec3<double>(pi_case::cos(tilt.phi) * pi_case::sin(tilt.theta),
                           pi_case::sin(tilt.phi) * pi_case::sin(tilt.theta),
                           pi_case::cos(tilt.theta));
}

vm::Vec3<int>
PliSimulator::LightDirectionComponent(const vm::Vec3<double> &dir_vec) const {
   vm::Vec3<int> dir_comp(0);
   for (auto itr = 0u; itr < 3; itr++) {
      if (dir_vec[itr] > 0)
         dir_comp[itr] = 1;
      else if (dir_vec[itr] < 0)
         dir_comp[itr] = -1;
   }
   return dir_comp;
}

std::vector<setup::Coordinates>
PliSimulator::CalcStartingLightPositions(const setup::Tilting &tilt) {
   // TODO: combine CalcStartingLightPositionsTilted and
   // CalcStartingLightPositionsUntilted
   if (setup_->untilt_sensor_view)
      return CalcStartingLightPositionsUntilted(tilt);
   else
      return CalcStartingLightPositionsTilted(tilt);
}

std::vector<setup::Coordinates>
PliSimulator::CalcStartingLightPositionsTilted(const setup::Tilting &tilt) {
   std::vector<setup::Coordinates> light_positions;
   /**
    * This function calculates the position of the initial beam on the bottom
    * plane. Layer of the fabric/label_field_.
    * 1. rotate the tissue in the LS (laboratory system) to the inclined
    * position.
    * 2. intersection point of the ccd position with the top plane
    * 3. Rotate coordinate to RF (Rotating Frame of tissue).
    * 4. light beam light beam is transferred to the lower tissue plane, the
    * starting point.
    */

   auto shift =
       LightDirectionUnitVector(tilt) / cos(tilt.theta) * dim_.global.z();

   // rotate tissue vector -> theta is inside tissue, refraction!
   // if light tiltes to one direction, the tissue tilts to the other
   const auto rot = vm::rot_pi_cases::Mat3RotZYZ(
       tilt.phi, asin(sin(-tilt.theta) * setup_->tissue_refrection), -tilt.phi);

   const auto rot_inv = vm::rot_pi_cases::Mat3RotZYZ(
       tilt.phi, -asin(sin(-tilt.theta) * setup_->tissue_refrection),
       -tilt.phi);

   // top plane: P = pc+pt+p1*s+p2*t
   // rotated top plane: P' = pc + rot(pt) + rot(p1) *s + rot(p2) * t
   const auto pc = vm::cast<double>(dim_.global) * 0.5;
   vm::Vec3<double> pt = {-pc.x(), -pc.y(), pc.z()};
   vm::Vec3<double> p1 = {1, 0, 0};
   vm::Vec3<double> p2 = {0, 1, 0};

   pt = vm::dot(rot, pt);
   p1 = vm::dot(rot, p1);
   p2 = vm::dot(rot, p2);

   // find position of light_beam(ccd_x=0, ccd_y=0) on top tissue plane
   // solve linear equation:
   const double a = p1.x();
   const double b = p2.x();
   const double c = p1.y();
   const double d = p2.y();

   // begin at ccd pos
   for (long long ccd_x = 0; ccd_x < dim_.global.x(); ccd_x++) {
      for (long long ccd_y = 0; ccd_y < dim_.global.y(); ccd_y++) {

         // light pos starts in the middle of pixel
         const double e = ccd_x + 0.5 - pc.x() - pt.x();
         const double f = ccd_y + 0.5 - pc.y() - pt.y();

         const double D = a * d - c * b;
         const double Ds = e * d - f * b;
         const double Dt = a * f - c * e;

         const double s = Ds / D;
         const double t = Dt / D;

         // TODO: test: elementwise inv-rot
         // position of ccd-tissue intersection in LS
         const vm::Vec3<double> p = pc + pt + p1 * s + p2 * t;

         // rotate back to RF
         const auto tissue = vm::dot(rot_inv, p - pc) + pc;

         // transform to bottom tissue position
         const double tis_x = tissue.x() - shift.x();
         const double tis_y = tissue.y() - shift.y();

         // check if processed by another mpi rank
         if (mpi_) {
            const auto coordinate = mpi_->coordinate();
            const auto global_coordinate = mpi_->global_coordinate();
            if (coordinate.x() != 0 && std::floor(tis_x) < dim_.offset.x())
               continue;
            if (coordinate.x() != global_coordinate.x() - 1 &&
                std::floor(tis_x) > (dim_.offset.x() + dim_.local.x() - 1))
               continue;
            if (coordinate.y() != 0 && std::floor(tis_y) < dim_.offset.y())
               continue;
            if (coordinate.y() != global_coordinate.y() - 1 &&
                std::floor(tis_y) > (dim_.offset.y() + dim_.local.y() - 1))
               continue;
         }

         setup::Coordinates data_point;
         data_point.tissue = {tis_x, tis_y, 0.0};
         data_point.ccd = {ccd_x, ccd_y};
         light_positions.push_back(data_point);
      }
   }

   if (light_positions.size() == 0) {
      if (mpi_)
         PyErr_WarnEx(PyExc_UserWarning,
                      ("rank " + std::to_string(mpi_->rank()) +
                       ": light_positions is empty")
                          .c_str(),
                      0);
      else
         PyErr_WarnEx(PyExc_UserWarning, "light_positions is empty", 0);
   }

   return light_positions;
}

std::vector<setup::Coordinates>
PliSimulator::CalcStartingLightPositionsUntilted(const setup::Tilting &tilt) {
   // TODO: check with outside tissue = background
   std::vector<setup::Coordinates> light_positions;

   auto shift =
       LightDirectionUnitVector(tilt) / cos(tilt.theta) * dim_.global.z();

   // begin at ccd pos
   for (long long ccd_x = 0; ccd_x < dim_.global.x(); ccd_x++) {
      for (long long ccd_y = 0; ccd_y < dim_.global.y(); ccd_y++) {

         // transform to bottom tissue position
         // light pos starts in the middle of pixel
         const double tis_x = (ccd_x + 0.5) - 0.5 * shift.x();
         const double tis_y = (ccd_y + 0.5) - 0.5 * shift.y();

         // check if processed by another process
         if (mpi_) {
            auto coordinate = mpi_->coordinate();
            auto global_coordinate = mpi_->global_coordinate();
            if (coordinate.x() != 0 && std::floor(tis_x) <= dim_.offset.x())
               continue;
            if (coordinate.x() != global_coordinate.x() - 1 &&
                std::floor(tis_x) > (dim_.offset.x() + dim_.local.x() - 1))
               continue;
            if (coordinate.y() != 0 && std::floor(tis_y) <= dim_.offset.y())
               continue;
            if (coordinate.y() != global_coordinate.y() - 1 &&
                std::floor(tis_y) > (dim_.offset.y() + dim_.local.y() - 1))
               continue;
         }

         setup::Coordinates data_point;
         data_point.tissue = {tis_x, tis_y, 0.0};
         data_point.ccd = {ccd_x, ccd_y};
         light_positions.push_back(data_point);
      }
   }

   if (light_positions.size() == 0) {
      if (mpi_)
         PyErr_WarnEx(PyExc_UserWarning,
                      ("rank " + std::to_string(mpi_->rank()) +
                       ": light_positions is empty")
                          .c_str(),
                      0);
      else
         PyErr_WarnEx(PyExc_UserWarning, "light_positions is empty", 0);
   }

   return light_positions;
}

// #############################################################################
// MPI functions
// #############################################################################

bool PliSimulator::CheckMPIHalo(const vm::Vec3<double> &local_pos,
                                const vm::Vec2<long long> &ccd_pos,
                                const vm::Vec3<int> &shift_direct,
                                const std::vector<vm::Vec4<double>> &s_vec) {

   if (!mpi_)
      return false;

   const auto low = dim_.offset;
   const auto coordinate = mpi_->coordinate();
   const auto global_coordinate = mpi_->global_coordinate();

   for (int i = 0; i < 3; i++) { // check x, y and z communication
      if ((local_pos[i] < 0.5 && shift_direct[i] == -1 && coordinate[i] > 0) ||
          (local_pos[i] > dim_.local[i] - 0.5 && shift_direct[i] == 1 &&
           coordinate[i] < global_coordinate[i] - 1)) {

         // ignore light in halo position at the beginning
         if (local_pos.z() == 0) {
            Abort(30000 + __LINE__); // should never happen
            // return true; // dont save, but do not process either
         }

#pragma omp critical
         mpi_->PushLightToBuffer(local_pos + vm::cast<double>(low), ccd_pos,
                                 s_vec, shift_direct[i] * (i + 1));
         return true;
      }
   }

   return false;
}

// #############################################################################
// only for testing
// #############################################################################

void PliSimulator::RunInterpolation(
    const vm::Vec3<long long> &dim, const vm::Vec3<long long> &dim_int,
    object::container::NpArray<int> label_field,
    object::container::NpArray<float> vector_field,
    object::container::NpArray<int> label_field_int,
    object::container::NpArray<float> vector_field_int,
    const setup::InterpMode interpolate) {

   dim_.local = dim;

   auto const scale = dim_int.x() / static_cast<double>(dim.x());

   if (dim.x() * dim.y() * dim.z() !=
       static_cast<long long>(label_field.size())) {
      std::cerr << "label_field.size(): " << label_field.size() << std::endl;
      Abort(30000 + __LINE__);
   }
   if (dim.x() * dim.y() * dim.z() * 3 !=
       static_cast<long long>(vector_field.size())) {
      std::cerr << "vector_field.size(): " << vector_field.size() << std::endl;
      Abort(30000 + __LINE__);
   }
   if (dim_int.x() * dim_int.y() * dim_int.z() !=
       static_cast<long long>(label_field_int.size())) {
      std::cerr << "label_field_int.size(): " << label_field_int.size()
                << std::endl;
      Abort(30000 + __LINE__);
   }
   if (dim_int.x() * dim_int.y() * dim_int.z() * 3 !=
       static_cast<long long>(vector_field_int.size())) {
      std::cerr << "vector_field_int.size(): " << vector_field_int.size()
                << std::endl;
      Abort(30000 + __LINE__);
   }
   // transfer data to class
   label_field_ =
       std::make_unique<object::container::NpArray<int>>(label_field);
   vector_field_ =
       std::make_unique<object::container::NpArray<float>>(vector_field);
   auto label_int =
       std::make_unique<object::container::NpArray<int>>(label_field_int);
   auto vector_int =
       std::make_unique<object::container::NpArray<float>>(vector_field_int);

   if (!setup_)
      setup_ = std::make_unique<setup::Setup>();

#pragma omp parallel for
   for (auto x = 0LL; x < dim_int.x(); x++) {
      for (auto y = 0LL; y < dim_int.y(); y++) {
         for (auto z = 0LL; z < dim_int.z(); z++) {
            auto const xs = (0.5 + x) / scale;
            auto const ys = (0.5 + y) / scale;
            auto const zs = (0.5 + z) / scale;

            const size_t idx =
                x * dim_int.y() * dim_int.z() + y * dim_int.z() + z;

            (*label_int)[idx] = GetLabel(llfloor(xs), llfloor(ys), llfloor(zs));

            auto const vec = GetVec(xs, ys, zs, interpolate);

            (*vector_int)[idx * 3 + 0] = vec[0];
            (*vector_int)[idx * 3 + 1] = vec[1];
            (*vector_int)[idx * 3 + 2] = vec[2];
         }
      }
   }
}
