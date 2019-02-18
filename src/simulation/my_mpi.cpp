#include "my_mpi.hpp"

#include <array>
#include <cassert>
#include <chrono>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <mpi.h>

#include "include/vemath.hpp"

MyMPI::MyMPI() {
   int flag = 0;
   MPI_Initialized(&flag);

   if (!flag) {
      ext_mpi_init_ = false;
      // std::cerr << "WARNING: calling MPI_Init(): " << flag << std::endl;
      // MPI_Init(NULL, NULL);
      // std::cerr << "WARNING: called MPI_Init()" << std::endl;
   }

   if (ext_mpi_init_) {
      MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_);
      MPI_Comm_size(MPI_COMM_WORLD, &numP_);
   } else {
      my_rank_ = 0;
      numP_ = 1;
   }

   if (debug_)
      std::cerr << "rank: " << my_rank_ << ", numP: " << numP_ << std::endl;

   // ClearBuffer();
}

MyMPI::~MyMPI() {
   // if (!ext_mpi_init_) {
   //    std::cout << "WARNING: calling MPI_Finalize()" << std::endl;
   //    MPI_Finalize();
   // }
}

void MyMPI::CreateCartGrid(const vm::Vec3<long long> global_dim) {

   const int max_dims = 3;
   const int reorder = 1;
   std::array<int, 3> period = {0, 0, 0};

   // calculate minimal communication volume
   vm::Vec3<long long> dim{};
   auto min_area = std::numeric_limits<long long>::max();

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
               global_coords_ = {x, y, z};
               min_area = area;
            }
         }
      }
   }

   if (ext_mpi_init_) {
      MPI_Cart_create(MPI_COMM_WORLD, max_dims, global_coords_.data(),
                      period.data(), reorder, &COMM_CART_);
      MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_);
      MPI_Cart_coords(COMM_CART_, my_rank_, 3, my_coords_.data());
   } else {
      global_coords_[0] = 1;
      global_coords_[0] = 1;
      global_coords_[0] = 1;
      my_coords_[0] = 0;
      my_coords_[1] = 0;
      my_coords_[2] = 0;
   }

   if (debug_)
      std::cout << "rank: " << my_rank_ << ", my_coords_:" << my_coords_
                << std::endl;

   CalcDimensions(global_dim);
}

void MyMPI::CalcDimensions(const vm::Vec3<long long> global_dim) {

   dim_vol_.global = global_dim;

   dim_vol_.local.x() =
       ceil(dim_vol_.global.x() / static_cast<double>(global_coords_[0]));
   dim_vol_.local.y() =
       ceil(dim_vol_.global.y() / static_cast<double>(global_coords_[1]));
   dim_vol_.local.z() =
       ceil(dim_vol_.global.z() / static_cast<double>(global_coords_[2]));

   vm::Vec3<long long> low{};
   vm::Vec3<long long> up{};

   low.x() = my_coords_.x() * dim_vol_.local.x();
   low.y() = my_coords_.y() * dim_vol_.local.y();
   low.z() = my_coords_.z() * dim_vol_.local.z();

   // +1 for halo
   up.x() = (my_coords_.x() + 1) * dim_vol_.local.x() - 1 + 1;
   up.y() = (my_coords_.y() + 1) * dim_vol_.local.y() - 1 + 1;
   up.z() = (my_coords_.z() + 1) * dim_vol_.local.z() - 1 + 1;

   assert(low.x() >= 0);
   assert(low.y() >= 0);
   assert(low.z() >= 0);

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
      std::cout << "rank: " << my_rank_ << ", my_coords_:" << my_coords_
                << std::endl;
      PrintDimensions(dim_vol_);
      std::cout << "rank " << my_rank_ << ": low:\t" << low << std::endl;
      std::cout << "rank " << my_rank_ << ": up:\t" << up << std::endl;
   }

   assert(dim_vol_.local.x() >= 0);
   assert(dim_vol_.local.y() >= 0);
   assert(dim_vol_.local.z() >= 0);

   // calculate ccd dimensions
   // dim_ccd_.low.x() = my_coords_.x() * dim_vol_.local.x();
   // dim_ccd_.up.x() = (my_coords_.x() + 1) * dim_vol_.local.x() - 1;
   // dim_ccd_.local.x() = dim_ccd_.up.x() - dim_ccd_.low.x() + 1;
   // dim_ccd_.global.x() = dim_vol_.global.x();

   // dim_ccd_.low.y() = my_coords_.y() * dim_vol_.local.y();
   // dim_ccd_.up.y() = (my_coords_.y() + 1) * dim_vol_.local.y() - 1;
   // dim_ccd_.local.y() = dim_ccd_.up.y() - dim_ccd_.low.y() + 1;
   // dim_ccd_.global.y() = dim_vol_.global.y();

   // dim_ccd_.low.z() = 0;
   // dim_ccd_.up.z() = 0;
   // dim_ccd_.local.z() = 0;
   // dim_ccd_.global.z() = 0;

   // if (debug)
   //    printDimensions(dim_ccd_);
}

void MyMPI::PrintDimensions(Dimensions dim) {
   std::cout << "rank " << my_rank_ << ": global:\t" << dim.global << std::endl;
   std::cout << "rank " << my_rank_ << ": local:\t" << dim.local << std::endl;
   std::cout << "rank " << my_rank_ << ": offset:\t" << dim.offset << std::endl;
   std::cout << "rank " << my_rank_ << ": origin:\t" << dim.origin << std::endl;
}

/*
void MyMPI::ClearBuffer() {

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

void MyMPI::PushLightToBuffer(vm::Vec3<double> pos, vm::Vec2<int> ccd,
                              std::vector<vm::Vec4<double>> light,
                              int direction) {

   if (direction == 0 || std::abs(direction) > 3) {
      CERROR << "Wrong direction value: " << direction << endl;
      exit(EXIT_FAILURE);
   }

   // index of sending direction
   int ind = std::distance(
       send_direction_.begin(),
       std::find(send_direction_.begin(), send_direction_.end(), direction));

   if (ind < 0 || ind > 5) {
      CERROR << "Wrong index: " << ind << ", direction: " << direction << endl;
      exit(EXIT_FAILURE);
   }

   // save all data in buffer. Order is important!
   for (auto elm : pos)
      snd_buffer_[ind].push_back(elm);

   for (auto elm : ccd)
      snd_buffer_[ind].push_back(static_cast<double>(elm));

   for (auto elm : light)
      for (auto e : elm)
         snd_buffer_[ind].push_back(e);

   if (num_rho_ != static_cast<int>(light.size())) {
      CERROR << "2*num_rho_ != light.size(): " << 2 * num_rho_
             << " != " << light.size() << endl;
      exit(EXIT_FAILURE);
   }
}

void MyMPI::CommunicateData(vector<Tissue2CCD> &scan_grid,
                            vector<vm::Vec4<double>> &intensity_buffer) {
   SndRcv();
   BufferToVariable(scan_grid, intensity_buffer);
   ClearBuffer();
}

void MyMPI::SndRcv() {
   // Start non blocking Isend:
   for (auto ind = 0u; ind < send_direction_.size(); ind++) {

      int direction = std::abs(send_direction_[ind]) - 1;
      int dispersion = send_direction_[ind] > 0 ? 1 : -1;
      int sender, receiver;

      MPI_Cart_shift(COMM_CART_, direction, dispersion, &sender, &receiver);
      if (receiver == MPI_PROC_NULL)
         continue;

      if (debug_)
         cout << "rank " << my_rank_ << " to " << receiver << " in direction "
              << direction << endl;

      MPI_Issend(snd_buffer_[ind].data(), snd_buffer_[ind].size(), MPI_DOUBLE,
                 receiver, tag_ + ind, COMM_CART_, &snd_rq_[ind]);
   }

   // Start blocking Recv:
   MPI_Status status;
   int count = 0;

   for (auto ind = 0u; ind < send_direction_.size(); ind++) {
      int direction = std::abs(send_direction_[ind]) - 1;
      int dispersion = send_direction_[ind] > 0 ? 1 : -1;
      int sender, receiver;

      MPI_Cart_shift(COMM_CART_, direction, dispersion, &sender, &receiver);
      if (sender == MPI_PROC_NULL)
         continue;

      MPI_Probe(sender, tag_ + ind, COMM_CART_, &status);
      MPI_Get_count(&status, MPI_DOUBLE, &count);
      if (count == MPI_UNDEFINED) {
         CERROR << "n == MPI_UNDEFINED" << endl;
         exit(EXIT_FAILURE);
      }
      rcv_buffer_[ind].resize(count);
      if (rcv_buffer_[ind].size() != static_cast<size_t>(count)) {
         CERROR << "resize didn't work" << endl;
         exit(EXIT_FAILURE);
      }
      MPI_Recv(rcv_buffer_[ind].data(), count, MPI_DOUBLE, MPI_ANY_SOURCE,
               tag_ + ind, COMM_CART_, &status);
   }
}

int MyMPI::Allreduce(int value) {
   int sum = 0;
   MPI_Allreduce(&value, &sum, 1, MPI_INT, MPI_SUM, COMM_CART_);
   return sum;
}

void MyMPI::BufferToVariable(vector<Tissue2CCD> &scan_grid,
                             vector<vm::Vec4<double>> &intensity_buffer) {
   scan_grid.clear();
   intensity_buffer.clear();

   if (num_rho_ <= 0) {
      CERROR << "num_rho_ <= 0" << endl;
      exit(EXIT_FAILURE);
   }

   Tissue2CCD pos_ccd;
   for (auto ind = 0u; ind < send_direction_.size(); ind++) {

      if (rcv_buffer_[ind].size() == 0)
         continue;

      for (auto i = 0u; i < rcv_buffer_[ind].size(); i += (5 + 4 * num_rho_)) {

         pos_ccd.tissue = {rcv_buffer_[ind][i], rcv_buffer_[ind][i + 1],
                           rcv_buffer_[ind][i + 2]};
         pos_ccd.ccd = {static_cast<int>(std::round(rcv_buffer_[ind][i + 3])),
                        static_cast<int>(std::round(rcv_buffer_[ind][i + 4]))};
         scan_grid.push_back(pos_ccd);

         size_t first_elm = 5;
         size_t last_elm = first_elm + 4 * num_rho_;

         for (auto j = first_elm; j < last_elm; j += 4) {
            intensity_buffer.push_back(vm::Vec4<double>(
                {{rcv_buffer_[ind][i + j + 0], rcv_buffer_[ind][i + j + 1],
                  rcv_buffer_[ind][i + j + 2], rcv_buffer_[ind][i + j + 3]}}));
         }
      }
   }
}
*/
