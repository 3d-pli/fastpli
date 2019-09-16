#include "simulator.hpp"

#include <exception>
#include <functional>
#include <iostream>
#include <memory>
#include <vector>

#include "helper.hpp"
#include "include/omp.hpp"
#include "include/vemath.hpp"
#include "my_mpi.hpp"
#include "objects/np_array_container.hpp"

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
   mpi_ = std::make_unique<MyMPI>(comm);
}

void PliSimulator::SetSetup(const setup::Setup setup) {
   setup_ = std::make_unique<setup::Setup>(setup);
}

void PliSimulator::CalculateDimensions(const vm::Vec3<long long> &global_dim) {

   if (global_dim.x() <= 0 || global_dim.y() <= 0 || global_dim.z() <= 0)
      throw std::invalid_argument("global_dim[any] <= 0: [" +
                                  std::to_string(global_dim.x()) + "," +
                                  std::to_string(global_dim.y()) + "," +
                                  std::to_string(global_dim.z()) + "]");

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
       label_field_.size())
      throw std::invalid_argument("dim_.local.x() * dim_.local.y() * "
                                  "dim_.local.z() != label_field_.size()");

   if (dim_.local.x() * dim_.local.y() * dim_.local.z() * 3ULL !=
       vector_field_.size())
      throw std::invalid_argument("dim_.local.x() * dim_.local.y() * "
                                  "dim_.local.z()* 3 != vector_field_.size()");
}

void PliSimulator::CheckDimension() {

   int min = std::numeric_limits<int>::max();
   int max = std::numeric_limits<int>::min();

   for (size_t i = 0; i < label_field_.size(); i++) {
      min = std::min(min, label_field_[i]);
      max = std::max(max, label_field_[i]);
   }

   if (min < 0 || max < 0)
      throw std::invalid_argument("label < 0 detected");

   if (static_cast<size_t>(max) >= properties_.size())
      throw std::invalid_argument("label exceed properties.size()");
}

std::vector<double>
PliSimulator::RunSimulation(const vm::Vec3<long long> &global_dim,
                            object::container::NpArray<int> label_field,
                            object::container::NpArray<float> vector_field,
                            setup::PhyProps properties, const double theta,
                            const double phi) {

   // transfer data to class
   label_field_ = std::move(label_field);
   vector_field_ = std::move(vector_field);
   properties_ = properties;

   // pre calculations and sanity checks
   if (!setup_)
      throw std::invalid_argument("pli_setup not set yet");
   CalculateDimensions(global_dim);
   CheckDimension();

   if (std::abs(theta) >= M_PI_2)
      throw std::invalid_argument("illegal light path");

   const auto n_rho = setup_->filter_rotations.size();
   const double lambda = setup_->wavelength * 1e-9;
   const double thickness = setup_->voxel_size * 1e-6 * setup_->step_size;

   // polarizer and lambda/4 retarder
   const double pol_x = 1; // TODO: via setup_
   const double pol_y = 1;

   vm::Mat4x4<double> polarizer_x(0);
   vm::Mat4x4<double> polarizer_y(0);
   vm::Mat4x4<double> m_lambda_4 = RetarderMatrix(M_PI_2, -M_PI_2);

   // FIXME: WRONG!! for p != 1
   polarizer_x(0, 0) = pol_x * pol_x;
   polarizer_x(0, 1) = pol_x * pol_x;
   polarizer_x(1, 0) = pol_x * pol_x;
   polarizer_x(1, 1) = pol_x * pol_x;

   polarizer_y(0, 0) = pol_y * pol_y;
   polarizer_y(0, 1) = -pol_y * pol_y;
   polarizer_y(1, 0) = -pol_y * pol_y;
   polarizer_y(1, 1) = pol_y * pol_y;

   polarizer_x *= 0.5;
   polarizer_y *= 0.5;

   // initial values
   vm::Vec4<double> signal_0 = {{setup_->light_intensity, 0, 0, 0}};
   signal_0 = vm::dot(m_lambda_4, vm::dot(polarizer_x, signal_0));
   const std::vector<vm::Vec4<double>> s_vec_0(n_rho, signal_0);

   const vm::Mat3x3<double> rotation =
       vm::rot_pi_cases::Mat3RotZYZ(phi, theta, -phi);

   const vm::Vec3<double> light_dir_vec = LightDirectionUnitVector(theta, phi);
   const vm::Vec3<int> light_dir_comp = LightDirectionComponent(light_dir_vec);
   const vm::Vec3<double> light_step = light_dir_vec * setup_->step_size;
   const auto TransformSensorPosToStart =
       GetSensorToStartTransformation(theta, phi);

   std::vector<double> intensity_signal(
       dim_.local.x() * dim_.local.y() * n_rho,
       std::numeric_limits<double>::quiet_NaN());

   auto light_positions = CalcStartingLightPositions(theta, phi);

   bool flag_all_done = false;
   while (!flag_all_done) {

#pragma omp parallel for
      for (size_t s = 0; s < light_positions.size(); s++) {

         auto light_position = light_positions[s];
         bool flag_ray_sended = false;

         auto local_pos = light_position.tissue - vm::cast<double>(dim_.offset);
         auto ccd_pos = light_position.ccd;
         auto s_vec = s_vec_0;

#pragma omp critical
         if (local_pos.z() >= 0.5) {
            s_vec.clear();
            for (auto i = 0u; i < n_rho; i++) {
               if (signal_buffer_.empty())
                  Abort(3111);

               s_vec.push_back(signal_buffer_.back());
               signal_buffer_.pop_back();
            }
            // reverse s_vec since signal_buffer_.pop_back();
            std::reverse(s_vec.begin(), s_vec.end());
         }

         for (; local_pos.z() > -0.5 && local_pos.z() < dim_.local.z() - 0.5;
              local_pos += light_step) {

            // check if communication is neccesary
            if (CheckMPIHalo(local_pos, light_dir_comp, s_vec, light_position))
               // next light ray
               break;

            auto label = GetLabel(local_pos);

            // calculate physical parameters
            auto dn = properties_.dn(label);
            auto mu = properties_.mu(label) * 1e3;
            const auto attenuation = pow(exp(-0.5 * mu * thickness), 2);

            if (dn == 0) {
               if (mu == 0)
                  continue;

               for (auto rho = 0u; rho < n_rho; rho++)
                  s_vec[rho] *= attenuation;
               continue;
            }

            auto vec = GetVec(local_pos, setup_->interpolate);
            const auto d_rel = dn * 4.0 * thickness / lambda;

            if (theta != 0)
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

               auto a1 = s_vec[rho][1] * cos_b - s_vec[rho][2] * sin_b;
               auto c1 = s_vec[rho][1] * sin_b + s_vec[rho][2] * cos_b;
               auto b1 = c1 * cos_r + s_vec[rho][3] * sin_r;

               // is equivalent to (R*M*R*S)*att
               s_vec[rho] = {{s_vec[rho][0], a1 * cos_b + b1 * sin_b,
                              -a1 * sin_b + b1 * cos_b,
                              -c1 * sin_r + s_vec[rho][3] * cos_r}};
               s_vec[rho] *= attenuation;
            }
         }

         if (!flag_ray_sended) {
            // save only, if light ray has reached the end of the volume
            size_t ccd_idx = (ccd_pos.x() - dim_.offset.x()) * dim_.local.y() +
                             (ccd_pos.y() - dim_.offset.y());

            // save signal only if light beam did not leave xy border
            if (local_pos.x() >= -0.5 && local_pos.y() >= -0.5 &&
                local_pos.x() < static_cast<double>(dim_.local.x()) - 0.5 &&
                local_pos.y() < static_cast<double>(dim_.local.y()) - 0.5) {

#ifndef NDEBUG
               if (ccd_idx * n_rho >= intensity_signal.size())
                  Abort(3112);
#endif
               for (auto rho = 0u; rho < n_rho; rho++) {
                  s_vec[rho] = vm::dot(polarizer_y, s_vec[rho]);
                  assert(!std::isnan(s_vec[rho][0]));
                  intensity_signal[ccd_idx * n_rho + rho] = s_vec[rho][0];
               }
            }
         }
      }

#pragma omp critical
      if (mpi_) {
#ifndef NDEBUG
         if (!signal_buffer_.empty())
            Abort(3113);
#endif
         mpi_->CommunicateData(light_positions, signal_buffer_);
         int num_communications =
             mpi_->Allreduce(light_positions.size() + signal_buffer_.size());

#ifndef NDEBUG
         if (signal_buffer_.size() != light_positions.size() * n_rho)
            Abort(3114);
#endif
         flag_all_done = num_communications == 0;
         if (debug_)
            std::cout << "rank " << mpi_->my_rank()
                      << ": num_communications: " << num_communications
                      << std::endl;
      } else {
         flag_all_done = true;
      }
   }

   return intensity_signal;
}

int PliSimulator::GetLabel(const long long x, const long long y,
                           const long long z) const {
   size_t idx = x * dim_.local.y() * dim_.local.z() + y * dim_.local.z() + z;
#ifndef NDEBUG
   if (idx >= label_field_.size())
      Abort(3115);
#endif
   return label_field_[idx];
}

vm::Vec3<double> PliSimulator::GetVec(const long long x, const long long y,
                                      const long long z) const {
   size_t idx =
       x * dim_.local.y() * dim_.local.z() * 3 + y * dim_.local.z() * 3 + z * 3;
#ifndef NDEBUG
   if (idx >= vector_field_.size())
      Abort(3116);
#endif
   return vm::Vec3<double>(vector_field_[idx], vector_field_[idx + 1],
                           vector_field_[idx + 2]);
}
vm::Vec3<double> PliSimulator::GetVec(const double x, const double y,
                                      const double z,
                                      const bool interpolate) const {
   return InterpolateVec(x, y, z, interpolate);
}

vm::Vec3<double> PliSimulator::InterpolateVec(const double x, const double y,
                                              const double z,
                                              const bool interpolate) const {

   // NN interpolate
   if (!interpolate)
      return GetVec(static_cast<long long>(std::round(x)),
                    static_cast<long long>(std::round(y)),
                    static_cast<long long>(std::round(z)));

   // Trilinear interpolate
   auto x0 = static_cast<long long>(std::floor(x));
   auto y0 = static_cast<long long>(std::floor(y));
   auto z0 = static_cast<long long>(std::floor(z));

   auto x1 = static_cast<long long>(std::ceil(x));
   auto y1 = static_cast<long long>(std::ceil(y));
   auto z1 = static_cast<long long>(std::ceil(z));

   if (x0 < 0)
      x0 = 0;
   if (y0 < 0)
      y0 = 0;
   if (z0 < 0)
      z0 = 0;

   if (x1 >= dim_.local.x())
      x1 = dim_.local.x() - 1;
   if (y1 >= dim_.local.y())
      y1 = dim_.local.y() - 1;
   if (z1 >= dim_.local.z())
      z1 = dim_.local.z() - 1;

   if (x0 == x1 && y0 == y1 && z0 == z1)
      return GetVec(x0, y0, z0);

   auto xd = (x - x0) / (x1 - x0);
   auto yd = (y - y0) / (y1 - y0);
   auto zd = (z - z0) / (z1 - z0);

   if (x0 == x1)
      xd = 0;
   if (y0 == y1)
      yd = 0;
   if (z0 == z1)
      zd = 0;

   auto c000 = GetVec(x0, y0, z0);
   auto c100 = GetVec(x1, y0, z0);
   auto c010 = GetVec(x0, y1, z0);
   auto c110 = GetVec(x1, y1, z0);
   auto c001 = GetVec(x0, y0, z1);
   auto c101 = GetVec(x1, y0, z1);
   auto c011 = GetVec(x0, y1, z1);
   auto c111 = GetVec(x1, y1, z1);

   auto c00 = c000 * (1 - xd) + c100 * xd;
   auto c01 = c001 * (1 - xd) + c101 * xd;
   auto c10 = c010 * (1 - xd) + c110 * xd;
   auto c11 = c011 * (1 - xd) + c111 * xd;

   auto c0 = c00 * (1 - yd) + c10 * yd;
   auto c1 = c01 * (1 - yd) + c11 * yd;

   return c0 * (1 - zd) + c1 * zd;
}

vm::Vec3<double>
PliSimulator::LightDirectionUnitVector(const double theta,
                                       const double phi) const {

   // special tilting angles
   if (theta == 0)
      return vm::Vec3<double>(0, 0, 1);
   else if (phi == 0)
      return vm::Vec3<double>(sin(theta), 0.0, cos(theta));
   else if (phi == M_PI_2)
      return vm::Vec3<double>(0.0, sin(theta), cos(theta));
   else if (phi == M_PI)
      return vm::Vec3<double>(-sin(theta), 0.0, cos(theta));
   else if (phi == 3 * M_PI_2)
      return vm::Vec3<double>(0.0, -sin(theta), cos(theta));

   return vm::Vec3<double>(cos(phi) * sin(theta), sin(phi) * sin(theta),
                           cos(theta));
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

std::function<vm::Vec3<double>(long long, long long)>
PliSimulator::GetSensorToStartTransformation(const double theta,
                                             const double phi) const {

   if (setup_->untilt_sensor) {
      const auto dz = dim_.local.z() - 1;
      // FIXME: const auto dz = dim_.z() - 0.5;

      // special tilting angles
      vm::Vec2<double> shift;
      if (theta == 0)
         shift = vm::Vec2<double>(0.0, 0.0);
      else if (phi == 0)
         shift = vm::Vec2<double>(tan(theta) * dz, 0.0);
      else if (phi == M_PI_2)
         shift = vm::Vec2<double>(0.0, tan(theta) * dz);
      else if (phi == M_PI)
         shift = vm::Vec2<double>(-tan(theta) * dz, 0.0);
      else if (phi == 3 * M_PI_2)
         shift = vm::Vec2<double>(0.0, -tan(theta) * dz);
      else
         shift =
             vm::Vec2<double>(tan(theta) * cos(phi), tan(theta) * sin(phi)) *
             dz;

      return [=](long long x, long long y) -> vm::Vec3<double> {
         return vm::Vec3<double>(x - shift.x(), y - shift.y(), 0.0);
      };
   } else {
      // TODO: implement back projection
      return std::function<vm::Vec3<double>(long long, long long)>();
   }
}

std::vector<Coordinates>
PliSimulator::CalcStartingLightPositions(const double phi, const double theta) {

   std::vector<Coordinates> light_positions;

   Coordinates data_point;
   vm::Vec3<double> shift =
       LightDirectionUnitVector(phi, theta) * (dim_.global.z() - 1);

   // begin at ccd pos
   for (long long ccd_x = 0; ccd_x < dim_.global.x(); ccd_x++) {
      for (long long ccd_y = 0; ccd_y < dim_.global.y(); ccd_y++) {
         // transform to bottom position
         auto tis_x = ccd_x - shift.x();
         auto tis_y = ccd_y - shift.y();

         if (tis_x < dim_.offset.x() - 0.5 ||
             tis_x > (dim_.offset.x() + dim_.local.x() - 0.5))
            continue;
         if (tis_y < dim_.offset.y() - 0.5 ||
             tis_y > (dim_.offset.y() + dim_.local.y() - 0.5))
            continue;

         data_point.tissue = {tis_x, tis_y, 0.0};
         data_point.ccd = {ccd_x, ccd_y};

         light_positions.push_back(data_point);
      }
   }

   if (light_positions.size() == 0) {
      if (mpi_)
         std::cout << mpi_->my_rank() << ": Warning, light_positions is empty"
                   << std::endl;
      else
         std::cout << 0 << ": Warning, light_positions is empty" << std::endl;
   }

   return light_positions;
}

vm::Mat4x4<double> PliSimulator::RetarderMatrix(const double beta,
                                                const double ret) const {
   vm::Mat4x4<double> M = {
       {1, 0, 0, 0, 0, pow(cos(beta), 2) + cos(ret) * pow(sin(beta), 2),
        (1 - cos(ret)) * sin(beta) * cos(beta), -sin(ret) * sin(beta), 0,
        (1 - cos(ret)) * sin(beta) * cos(beta),
        pow(sin(beta), 2) + cos(ret) * pow(cos(beta), 2), sin(ret) * cos(beta),
        0, sin(ret) * sin(beta), -sin(ret) * cos(beta), cos(ret)}};
   return M;
}

bool PliSimulator::CheckMPIHalo(const vm::Vec3<double> &local_pos,
                                const vm::Vec3<int> &shift_direct,
                                const std::vector<vm::Vec4<double>> &s_vec,
                                const Coordinates &startpos) {

   if (!mpi_)
      return false;

   // this can happen, because CalcStartingLightPositions() is not sensitiv to
   // the mpi-halo yet
   if (local_pos.z() == 0)
      return true;

   const auto low = dim_.offset;
   const auto up = dim_.offset + dim_.local;

   bool flag = false;
   // x - halo communication:
   if ((local_pos.x() < 0 && shift_direct.x() == -1 && low.x() != 0) ||
       (local_pos.x() > dim_.local.x() - 0.5 && shift_direct.x() == 1 &&
        up.x() != dim_.global.x() - 0.5)) {
#pragma omp critical
      mpi_->PushLightToBuffer(local_pos + vm::cast<double>(low), startpos.ccd,
                              s_vec, shift_direct.x() * 1);
      flag = true;

      // y-halo communication:
   } else if ((local_pos.y() < 0 && shift_direct.y() == -1 && low.y() != 0) ||
              (local_pos.y() > dim_.local.y() - 0.5 && shift_direct.y() == 1 &&
               up.y() != dim_.global.y() - 0.5)) {
#pragma omp critical
      mpi_->PushLightToBuffer(local_pos + vm::cast<double>(low), startpos.ccd,
                              s_vec, shift_direct.y() * 2);
      flag = true;

      // z-halo communication:
   } else if ((local_pos.z() < 0 && shift_direct.z() == -1 && low.z() != 0) ||
              (local_pos.z() > dim_.local.z() - 0.5 && shift_direct.z() == 1 &&
               up.z() != dim_.global.z() - 0.5)) {
#pragma omp critical
      mpi_->PushLightToBuffer(local_pos + vm::cast<double>(low), startpos.ccd,
                              s_vec, shift_direct.z() * 3);
      flag = true;
   }

   return flag;
}
