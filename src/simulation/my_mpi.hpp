#ifndef SIMULATION_MY_MPI_HPP_
#define SIMULATION_MY_MPI_HPP_

#include <array>
#include <tuple>
#include <vector>

#include "include/vemath.hpp"
#include "setup.hpp"

#include <mpi.h>

class MyMPI {
 public:
   MyMPI(const MPI_Comm comm);
   ~MyMPI();

   void CreateCartGrid(vm::Vec3<long long> dim_global);

   // communication
   void set_num_rho(int num_rho) { num_rho_ = num_rho; };

   void PushLightToBuffer(vm::Vec3<double> pos, vm::Vec2<long long> ccd,
                          std::vector<vm::Vec4<double>> light, int direction);

   std::tuple<std::vector<setup::Coordinates>, std::vector<vm::Vec4<double>>>
   CommunicateData();

   int Allreduce(int value);

   // getter
   MPI_Comm comm() const { return my_comm_; };
   int my_rank() const { return my_rank_; };
   int num_of_proz() const { return numP_; };
   vm::Vec3<int> my_coords() const { return my_coords_; };
   vm::Vec3<int> global_coords() const { return global_coords_; };

   setup::Dimensions dim_vol() const { return dim_vol_; };

 private:
   void SndRcv();
   void ClearInternalBuffer();
   std::tuple<std::vector<setup::Coordinates>, std::vector<vm::Vec4<double>>>
   InternalBufferToVariable();

#ifndef NDEBUG
   const bool debug_ = true;
#else
   const bool debug_ = false;
#endif

   MPI_Comm my_comm_;
   int my_rank_ = 0;
   int numP_ = 1;

   setup::Dimensions dim_vol_{};

   vm::Vec3<int> my_coords_ = {-1, -1, -1};
   vm::Vec3<int> global_coords_ = {-1, -1, -1};
   MPI_Comm COMM_CART_;

   int num_rho_ = 0;
   const int tag_ = 0;
   const std::array<int, 6> send_direction_{{1, 2, 3, -1, -2, -3}};

   int number_of_dests_ = 0;
   int number_of_recvs_ = 0;
   std::array<int, 6> rcv_flag_;       // 6 directions
   std::array<MPI_Request, 6> snd_rq_; // initialization in constructor
   std::array<std::vector<double>, 6> rcv_buffer_;
   std::array<std::vector<double>, 6> snd_buffer_;

   void CalcDimensions(const vm::Vec3<long long> dim);
};

#endif // SIMULATION_MY_MPI_HPP_
