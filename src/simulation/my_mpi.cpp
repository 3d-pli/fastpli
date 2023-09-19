#include "my_mpi.hpp"

#include <array>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <vector>

#include <mpi.h>

#include "include/vemath.hpp"

MyMPI::MyMPI(const MPI_Comm comm) {
   my_comm_ = comm;
   MPI_Comm_rank(my_comm_, &rank_);
   MPI_Comm_size(my_comm_, &numP_);

   if (debug_)
      std::cout << "rank: " << rank_ << ", numP: " << numP_
                << ", comm: " << comm << std::endl;
}

MyMPI::~MyMPI() {}

void MyMPI::CreateCartGrid(const vm::Vec3<int64_t> global_dim) {

   const int max_dims = 3;
   const int reorder = 1;
   std::array<int, 3> period = {0, 0, 0};

   // calculate minimal communication volume
   vm::Vec3<int64_t> dim{};
   auto min_area = std::numeric_limits<int64_t>::max();

   // for(auto z = 1; z <= this->numP; z++){
   for (auto z = 1; z <= 1; z++) { // z dominates communication
      for (auto y = 1; y <= numP_; y++) {
         for (auto x = 1; x <= numP_; x++) {
            if (x * y * z != numP_)
               continue;

            dim.x() = ceil(global_dim.x() / static_cast<double>(x));
            dim.y() = ceil(global_dim.y() / static_cast<double>(y));
            dim.z() = ceil(global_dim.z() / static_cast<double>(z));

            // "area" for no z communication
            auto area = dim.x() + dim.y();

            if (area < min_area) {
               global_coordinate_ = {x, y, z};
               min_area = area;
            }
         }
      }
   }

   MPI_Cart_create(my_comm_, max_dims, global_coordinate_.data(), period.data(),
                   reorder, &COMM_CART_);
   MPI_Comm_rank(my_comm_, &rank_);
   MPI_Cart_coords(COMM_CART_, rank_, 3, coordinate_.data());

   CalcDimensions(global_dim);
}

void MyMPI::CalcDimensions(const vm::Vec3<int64_t> global_dim) {

   dim_vol_.global = global_dim;

   dim_vol_.local.x() =
       ceil(dim_vol_.global.x() / static_cast<double>(global_coordinate_[0]));
   dim_vol_.local.y() =
       ceil(dim_vol_.global.y() / static_cast<double>(global_coordinate_[1]));
   dim_vol_.local.z() =
       ceil(dim_vol_.global.z() / static_cast<double>(global_coordinate_[2]));

   vm::Vec3<int64_t> low{};
   vm::Vec3<int64_t> up{};

   low.x() = coordinate_.x() * dim_vol_.local.x();
   low.y() = coordinate_.y() * dim_vol_.local.y();
   low.z() = coordinate_.z() * dim_vol_.local.z();

   // +1 for halo
   up.x() = (coordinate_.x() + 1) * dim_vol_.local.x() - 1 + 1;
   up.y() = (coordinate_.y() + 1) * dim_vol_.local.y() - 1 + 1;
   up.z() = (coordinate_.z() + 1) * dim_vol_.local.z() - 1 + 1;

   if (vm::any_of(low, [](int64_t i) { return i < 0; })) {
      MPI_Abort(my_comm_, 10000 + __LINE__);
   }

   dim_vol_.offset = low;

   if (up.x() >= dim_vol_.global.x())
      up.x() = dim_vol_.global.x() - 1;
   if (up.y() >= dim_vol_.global.y())
      up.y() = dim_vol_.global.y() - 1;
   if (up.z() >= dim_vol_.global.z())
      up.z() = dim_vol_.global.z() - 1;
   dim_vol_.local = up - low + 1;

   // if there are to many threads, local dimension can get negativ
   // setting them to 0 so they are idle
   if (dim_vol_.local.x() < 0)
      dim_vol_.local.x() = 0;
   if (dim_vol_.local.y() < 0)
      dim_vol_.local.y() = 0;
   if (dim_vol_.local.z() < 0)
      dim_vol_.local.z() = 0;

   if (debug_) {
      std::cout << "rank " << rank_ << ": coordinate: " << coordinate_
                << std::endl;
      std::cout << "rank " << rank_ << ": low: " << low << std::endl;
      std::cout << "rank " << rank_ << ": up: " << up << std::endl;
   }

   if (vm::any_of(dim_vol_.local, [&](int64_t i) { return i < 0; })) {
      MPI_Abort(my_comm_, 10000 + __LINE__);
   }
}

void MyMPI::ClearInternalBuffer() {

   number_of_dests_ = 0;
   number_of_recvs_ = 0;

   for (auto &elm : snd_buffer_)
      elm.clear();
   for (auto &elm : rcv_buffer_)
      elm.clear();
   for (auto &elm : snd_rq_)
      elm = MPI_REQUEST_NULL;
   for (auto &elm : rcv_flag_)
      elm = 0;
}

void MyMPI::PushLightToBuffer(vm::Vec3<double> pos, vm::Vec2<int64_t> ccd,
                              std::vector<vm::Vec4<double>> light,
                              int direction) {

#ifndef NDEBUG
   if (direction == 0 || std::abs(direction) > 3)
      MPI_Abort(my_comm_, 10000 + __LINE__);
#endif

   // index of sending direction
   int ind = std::distance(
       send_direction_.begin(),
       std::find(send_direction_.begin(), send_direction_.end(), direction));

#ifndef NDEBUG
   if (ind < 0 || ind >= 6)
      MPI_Abort(my_comm_, 10000 + __LINE__);
#endif

   // save all data in buffer. Order is important!
   for (auto &elm : pos)
      snd_buffer_[ind].push_back(elm);

   for (auto &elm : ccd)
      snd_buffer_[ind].push_back(*reinterpret_cast<double *>(&elm));

   for (auto elm : light)
      for (auto e : elm)
         snd_buffer_[ind].push_back(e);

#ifndef NDEBUG
   if (num_rho_ != static_cast<int>(light.size()))
      MPI_Abort(my_comm_, 10000 + __LINE__);
#endif
}

std::tuple<std::vector<setup::Coordinates>, std::vector<vm::Vec4<double>>>
MyMPI::CommunicateData() {
   SndRcv();
   auto data = InternalBufferToVariable();
   ClearInternalBuffer();
   return data;
}

void MyMPI::SndRcv() {
   int sender, receiver;

   // Start non blocking Isend:
   for (auto ind = 0u; ind < send_direction_.size(); ind++) {

      int direction = std::abs(send_direction_[ind]) - 1;
      if (send_direction_[ind] == 0)
         continue;

      int displacement = send_direction_[ind] > 0 ? 1 : -1;

      MPI_Cart_shift(COMM_CART_, direction, displacement, &sender, &receiver);
      if (receiver == MPI_PROC_NULL)
         continue;

      if (debug_)
         std::cout << "rank " << rank_ << " to " << receiver << " in direction "
                   << direction << ": size: " << snd_buffer_[ind].size()
                   << std::endl;

      MPI_Issend(snd_buffer_[ind].data(), snd_buffer_[ind].size(), MPI_DOUBLE,
                 receiver, tag_ + ind, COMM_CART_, &snd_rq_[ind]);
   }

   // Start blocking Recv:
   MPI_Status status;
   int count = 0;

   for (auto ind = 0u; ind < send_direction_.size(); ind++) {
      int direction = std::abs(send_direction_[ind]) - 1;
      if (send_direction_[ind] == 0)
         continue;
      int displacement = send_direction_[ind] > 0 ? 1 : -1;

      MPI_Cart_shift(COMM_CART_, direction, displacement, &sender, &receiver);
      if (sender == MPI_PROC_NULL)
         continue;

      MPI_Probe(sender, tag_ + ind, COMM_CART_, &status);
      MPI_Get_count(&status, MPI_DOUBLE, &count);
      if (count == MPI_UNDEFINED)
         MPI_Abort(my_comm_, 10000 + __LINE__);

      if (debug_)
         std::cout << "rank " << rank_ << " from " << sender << " in direction "
                   << direction << ": size: " << count << std::endl;

      rcv_buffer_[ind].resize(count);
      MPI_Recv(rcv_buffer_[ind].data(), count, MPI_DOUBLE, MPI_ANY_SOURCE,
               tag_ + ind, COMM_CART_, &status);
   }
}

int MyMPI::AllreduceSum(int value) {
   int sum = 0;
   MPI_Allreduce(&value, &sum, 1, MPI_INT, MPI_SUM, COMM_CART_);
   return sum;
}

std::vector<int> MyMPI::AllreduceSum(std::vector<int> v) {
   std::vector<int> a(v.size());
   MPI_Allreduce(v.data(), a.data(), v.size(), MPI_INT, MPI_SUM, COMM_CART_);
   return a;
}

std::vector<double> MyMPI::AllreduceMax(std::vector<double> v) {
   std::vector<double> a(v.size());
   MPI_Allreduce(v.data(), a.data(), v.size(), MPI_DOUBLE, MPI_MAX, COMM_CART_);
   return a;
}

std::tuple<std::vector<setup::Coordinates>, std::vector<vm::Vec4<double>>>
MyMPI::InternalBufferToVariable() {
   std::vector<setup::Coordinates> light_positions;
   std::vector<vm::Vec4<double>> light_signals;

#ifndef NDEBUG
   if (num_rho_ <= 0)
      MPI_Abort(my_comm_, 10000 + __LINE__);
#endif

   setup::Coordinates pos_ccd;
   for (auto ind = 0u; ind < send_direction_.size(); ind++) {

      if (rcv_buffer_[ind].size() == 0)
         continue;

      for (auto i = 0u; i < rcv_buffer_[ind].size(); i += (5 + 4 * num_rho_)) {

         pos_ccd.tissue = {rcv_buffer_[ind][i], rcv_buffer_[ind][i + 1],
                           rcv_buffer_[ind][i + 2]};

         pos_ccd.ccd.x() =
             *reinterpret_cast<int64_t *>(&rcv_buffer_[ind][i + 3]);
         pos_ccd.ccd.y() =
             *reinterpret_cast<int64_t *>(&rcv_buffer_[ind][i + 4]);
         light_positions.push_back(pos_ccd);

         size_t first_elm = 5;
         size_t last_elm = first_elm + 4 * num_rho_;

         for (auto j = first_elm; j < last_elm; j += 4) {
            light_signals.push_back(vm::Vec4<double>(
                rcv_buffer_[ind][i + j + 0], rcv_buffer_[ind][i + j + 1],
                rcv_buffer_[ind][i + j + 2], rcv_buffer_[ind][i + j + 3]));
         }
      }
#ifndef NDEBUG
      if (light_positions.size() * num_rho_ != light_signals.size())
         MPI_Abort(my_comm_, 10000 + __LINE__);
#endif
   }

   return std::make_tuple(light_positions, light_signals);
}
