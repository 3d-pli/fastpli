#ifndef SIMULATION_MY_MPI_HPP_
#define SIMULATION_MY_MPI_HPP_

#include <array>
#include <vector>

#include "helper.hpp"
#include "include/vemath.hpp"

#include <mpi.h>

class MyMPI {
 public:
   MyMPI();
   ~MyMPI();

   // MPI functions
   void Barrier() { MPI_Barrier(MPI_COMM_WORLD); };
   void CreateCartGrid(vm::Vec3<long long> dim_global);

   // TODO: mpi for simulation.cpp
   // communication
   // void set_num_rho(int num_rho) { num_rho_ = num_rho; };
   // void ClearBuffer();
   // void PushLightToBuffer(vm::Vec3<double> pos, vm::Vec2<int> ccd,
   //                        std::vector<vm::Vec4<double>> light, int
   //                        direction);

   // void CommunicateData(std::vector<Tissue2CCD> &scan_xy,
   //                      std::vector<vm::Vec4<double>> &intensity_buffer);
   // void SndRcv();
   // void BufferToVariable(std::vector<Tissue2CCD> &scan_xy,
   //                       std::vector<vm::Vec4<double>> &intensity_buffer);
   // int Allreduce(int value);

   // getter
   int my_rank() const { return my_rank_; };
   vm::Vec3<int> my_coords() const { return my_coords_; };
   vm::Vec3<int> global_coords() const { return global_coords_; };
   int num_of_proz() const { return numP_; };

   Dimensions dim_vol() const { return dim_vol_; };
   Dimensions dim_ccd() const { return dim_ccd_; };
   void PrintDimensions(Dimensions dim);

 private:
   bool ext_mpi_init_ = true;
   int my_rank_ = -1;
   int numP_ = 0;
   const bool debug_ = false;

   // int num_rho_ = 0;

   Dimensions dim_vol_{};
   Dimensions dim_ccd_{};

   vm::Vec3<int> my_coords_ = {-1, -1, -1};
   vm::Vec3<int> global_coords_ = {-1, -1, -1};
   MPI_Comm COMM_CART_;

   // const int tag_ = 0;
   // const std::array<int, 6> send_direction_{{1, 2, 3, -1, -2, -3}};

   // int number_of_dests_ = 0;
   // int number_of_recvs_ = 0;
   // std::array<int, 6> rcv_flag_;       // 6 directions
   // std::array<MPI_Request, 6> snd_rq_; // initialization in constructor
   // std::array<std::vector<double>, 6> rcv_buffer_;
   // std::array<std::vector<double>, 6> snd_buffer_;

   void CalcDimensions(const vm::Vec3<long long> dim);
};

#endif // SIMULATION_MY_MPI_HPP_
