/// Class representing the problem to be solved. It contains the grid, the time step, the final time, the initial condition, the boundary conditions, the source term, the exact solution, the numerical solution, the error, and the output file.
/// the boundary conditions, the source term, the exact solution, the numerical solver.

#ifndef CORE_HPP
#define CORE_HPP

#include "utils.hpp"
#include "boundary.hpp"
#include "grid.hpp"
#include <string>
#include <cmath>
#include <filesystem>
#include <mpi.h>

class IcoNS
{
public:

  IcoNS(const std::string &input_file, const std::string &output_file)
      : input_file(input_file),
        output_file(output_file)

 
  {
  }

  void preprocessing(/*std::string &input_file*/); // grid initialization.


  Real functionF_u(const std::array<Real, NX *(NY + 1) * (NZ + 1)> &u, const std::array<Real, (NX + 1) * NY *(NZ + 1)> &v, const std::array<Real, (NX + 1) * (NY + 1) * NZ> &w, size_t i, size_t j, size_t k, Real t); // compute the source term.
  Real functionF_v(const std::array<Real, NX *(NY + 1) * (NZ + 1)> &u, const std::array<Real, (NX + 1) * NY *(NZ + 1)> &v, const std::array<Real, (NX + 1) * (NY + 1) * NZ> &w, size_t i, size_t j, size_t k, Real t); // compute the source term.
  Real functionF_w(const std::array<Real, NX *(NY + 1) * (NZ + 1)> &u, const std::array<Real, (NX + 1) * NY *(NZ + 1)> &v, const std::array<Real, (NX + 1) * (NY + 1) * NZ> &w, size_t i, size_t j, size_t k, Real t); // compute the source term.
  Real functionG_u(size_t i, size_t j, size_t k, Real t);                                                                                                                                                                    // compute the source term.
  Real functionG_v(size_t i, size_t j, size_t k, Real t);                                                                                                                                                                    // compute the source term.
  Real functionG_w(size_t i, size_t j, size_t k, Real t);                                                                                                                                                                    // compute the source term.


  void solve();                      // solve the problem saving the ouput.
  void solve_time_step(Real time); // solve a time step.

  Real error_comp_X(const Real t);
  Real error_comp_Y(const Real t);
  Real error_comp_Z(const Real t);
  Real L2_error(const Real t); // compute the L2 norm

  void output_u(const std::string& filename, Grid& print); // write the output file.
  void output_v(const std::string& filename, Grid& print); // write the output file.
  void output_w(const std::string& filename, Grid& print); // write the output file.

  void setParallelization();
  void IcoNS::exchangeData(std::array& grid_loc, int& newDimX, int newDimY,int dim_z, MPI_Datatype MPI_face_x, MPI_Datatype MPI_face_y);
private:
  Grid grid; // grid of the domain.
  Boundary boundary;
  ExactSolution exact_solution;
  std::array<Real, (NX * (NY + 1) * (NZ + 1))> Y2_x{};
  std::array<Real, ((NX + 1) * NY * (NZ + 1))> Y2_y{};
  std::array<Real, ((NX + 1) * (NY + 1) * NZ)> Y2_z{};
  std::array<Real, (NX * (NY + 1) * (NZ + 1))> Y3_x{};
  std::array<Real, ((NX + 1) * NY * (NZ + 1))> Y3_y{};
  std::array<Real, ((NX + 1) * (NY + 1) * NZ)> Y3_z{};
  std::string input_file;  // input file.
  std::string output_file; // output file.

  //Parallel vars
  int rank, size;
  MPI_Comm cart_comm;
  int dims[2] = {PX, PY};
  int periods[2] = {1,1};

  //Left boundary for x dim
  //since i don't know if i have to count boundary for each process, and if to count them from left or right i'll use these vars as booleans
  int lbx=0, rbx=0, lby =0, rby=0;


  int coords[2];
  int neighbors[4];

  std::array<Real,newDimX_x*newDimY_x*dim_z> grid_loc_x{};
  std::array<Real,newDimX_y*newDimY_y*dim_z> grid_loc_y{};
  std::array<Real,newDimX_z*newDimY_z*dim_z> grid_loc_z{};

  MPI_Datatype MPI_face_x_x;
  MPI_Datatype MPI_face_x_y;
  MPI_Datatype MPI_face_x_z;

  MPI_Datatype MPI_face_y_x;
  MPI_Datatype MPI_face_y_y;
  MPI_Datatype MPI_face_y_z;

  MPI_Request req1;
  MPI_Request req2;
  MPI_Request req3;
  MPI_Request req4;

};

#endif // CORE_HPP
