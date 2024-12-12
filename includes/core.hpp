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
  IcoNS(MPI_Comm cart_comm, const std::string &input_file, const std::string &output_file, int rank, int size)
      : cart_comm(cart_comm),
        input_file(input_file),
        output_file(output_file),
        rank(rank),
        size(size), 
        dim_x_x(NX / PX),
        dim_y_x((NY + 1) / PY),
        dim_x_y((NX + 1) / PX),
        dim_y_y(NY / PY),
        dim_x_z((NX + 1) / PX),
        dim_y_z((NY + 1) / PY),
        dim_z(NZ + 1),
        dim_z_z(NZ)
  {
    other_dim_x_x = dim_x_x;
    other_dim_y_x = dim_y_x;
    other_dim_x_y = dim_x_y;
    other_dim_y_y = dim_y_y;
    other_dim_x_z = dim_x_z;
    other_dim_y_z = dim_y_z;

    // grid_loc_x[newDimX_x * newDimY_x * (NZ + 1), 0.0];
    // grid_loc_y[newDimX_y * newDimY_y * (NZ + 1), 0.0];
    // grid_loc_z[newDimX_z * newDimY_z * (NZ), 0.0];

    // Y2_x[newDimX_x * newDimY_x * (NZ + 1), 0.0];
    // Y2_y[newDimX_y * newDimY_y * (NZ + 1), 0.0];
    // Y2_z[newDimX_z * newDimY_z * (NZ), 0.0];

    // Y3_x[newDimX_x * newDimY_x * (NZ + 1), 0.0];
    // Y3_y[newDimX_y * newDimY_y * (NZ + 1), 0.0];
    // Y3_z[newDimX_z * newDimY_z * (NZ), 0.0];
  }
  /**
  * @brief The method is called by the program to initialize the grid.
  *
  * @param input_files 
  */
  void preprocessing(/*std::string &input_file*/);

  /**
  * @brief The method is called by the program to initialize the parallelization paradigm.
  *
  * @param None
  */
  void setParallelization();

  /**
  * @brief The method is called by the program to exchange the ghost points between the processes.
  *
  * @param grid_loc 
  * @param newDimX 
  * @param grid_loc 
  */
  void exchangeData(std::vector<Real> &grid_loc, int newDimX, int newDimY, int dim_z, MPI_Datatype MPI_face_x, MPI_Datatype MPI_face_y);

  /**
  * @brief The method is called by the program multiple during the time step, 
  *        in order to update the values of the boundaries at each
  *        requested time t, calculating the approximated ones too.
  *
  * @param 
    */
  Real functionF_u(const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &w, int i, int j, int k, Real t); 
  Real functionF_v(const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &w, int i, int j, int k, Real t); 
  Real functionF_w(const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &w, int i, int j, int k, Real t); 
  
  // Manufactered solution functions.
  Real functionG_u(int i, int j, int k, Real t);                                                                                     
  Real functionG_v(int i, int j, int k, Real t);                                                                                     
  Real functionG_w(int i, int j, int k, Real t);                                                                                     

  void solve();                    
  void solve_time_step(Real time);

  // L2 error functions.
  Real error_comp_X(const Real t);
  Real error_comp_Y(const Real t);
  Real error_comp_Z(const Real t);
  Real L2_error(const Real t); 

  // TODO write the output file.

  

private:

  // Parallel variables
  int rank, size;
  MPI_Comm cart_comm;
  int dims[2] = {PX, PY};
  int periods[2] = {1, 1};

  Boundary boundary;
  ExactSolution exact_solution;
  
  std::vector<Real> Y2_x{};
  std::vector<Real> Y2_y{};
  std::vector<Real> Y2_z{};
  std::vector<Real> Y3_x{};
  std::vector<Real> Y3_y{};
  std::vector<Real> Y3_z{};
  
  std::string input_file;  // input file.
  std::string output_file; // output file.

  // Left boundary for x dim
  // since i don't know if i have to count boundary for each process, and if to count them from left or right i'll use these vars as booleans
  int lbx = 0, rbx = 0, lby = 0, rby = 0;

  int coords[2];
  int neighbors[4];

  std::vector<Real> grid_loc_x{};
  std::vector<Real> grid_loc_y{};
  std::vector<Real> grid_loc_z{};

  MPI_Datatype MPI_face_x_x, MPI_face_y_x;
  MPI_Datatype MPI_face_x_y, MPI_face_y_y;
  MPI_Datatype MPI_face_x_z, MPI_face_y_z;

  MPI_Request reqs[4];

  int dim_x_x, dim_y_x;
  int dim_x_y, dim_y_y, dim_y_z;
  int dim_z, dim_x_z, dim_z_z;

  int newDimX_x, newDimY_x; 
  int newDimX_y, newDimY_y;
  int newDimX_z, newDimY_z;

  int other_dim_x_x, other_dim_y_x;
  int other_dim_x_y, other_dim_y_y;
  int other_dim_x_z, other_dim_y_z;
};

#endif // CORE_HPP
