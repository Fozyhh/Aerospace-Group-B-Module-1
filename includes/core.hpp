/**
 * @file core.hpp
 * @brief Core implementation of the incompressible Navier-Stokes solver
 *
 * This file contains the main solver class for incompressible Navier-Stokes equations.
 * It handles grid management, parallel computation, boundary conditions, and time-stepping
 * solutions using MPI for parallel processing.
 */

#ifndef CORE_HPP
#define CORE_HPP

#define Real double

#include "utils.hpp"
#include "boundary.hpp"
#include "grid.hpp"
#include <string>
#include <cmath>
#include <filesystem>
#include <mpi.h>

/**
 * @class IcoNS
 * @brief Incompressible Navier-Stokes solver class
 *
 * Handles the numerical solution of incompressible Navier-Stokes equations
 * including parallel domain decomposition, boundary conditions, and time evolution.
 */
class IcoNS
{
public:
  /**
   * @brief Constructor for the IcoNS solver
   * @param cart_comm MPI Cartesian communicator
   * @param input_file Path to input configuration file
   * @param output_file Path to output results file
   * @param rank MPI process rank
   * @param size Total number of MPI processes
   */
  IcoNS(MPI_Comm cart_comm, const std::string &input_file, const std::string &output_file, int rank, int size)
      : cart_comm(cart_comm),
        input_file(input_file),
        output_file(output_file),
        rank(rank),
        size(size)
  {
    parse_input(input_file);

    dims[0] = PX;
    dims[1] = PY;

    dim_x_x = NX / PX;
    dim_y_x = (NY + 1) / PY;
    dim_x_y = (NX + 1) / PX;
    dim_y_y = NY / PY;
    dim_x_z = (NX + 1) / PX;
    dim_y_z = (NY + 1) / PY;
    dim_z = NZ + 1;
    dim_z_z = NZ;

    other_dim_x_x = dim_x_x;
    other_dim_y_x = dim_y_x;
    other_dim_x_y = dim_x_y;
    other_dim_y_y = dim_y_y;
    other_dim_x_z = dim_x_z;
    other_dim_y_z = dim_y_z;
  }

  /**
   * @brief Initializes the grid and problem setup
   */
  void preprocessing(/*std::string &input_file*/);

  /**
   * @brief Sets up parallel communication patterns and domain decomposition
   */
  void setParallelization();

  /**
   * @brief Exchanges ghost cell data between neighboring processes
   * @param grid_loc Local grid data
   * @param newDimX New X dimension after decomposition
   * @param newDimY New Y dimension after decomposition
   * @param dim_z Z dimension
   * @param MPI_face_x MPI datatype for X-direction face exchange
   * @param MPI_face_y MPI datatype for Y-direction face exchange
   * @param sameX Only for periodic boundary, skip identical datas during communication(set to 0 only for communication x grids)
   * @param sameY Only for periodic boundary, skip identical datas during communication(set to 0 only for communication y grids)
   */
  void exchangeData(std::vector<Real> &grid_loc, int newDimX, int newDimY, int dim_z, MPI_Datatype MPI_face_x, MPI_Datatype MPI_face_y,int sameX, int sameY);

  /**
   * @brief Computes the F function for u-velocity component
   * @param u X-velocity vector
   * @param v Y-velocity vector
   * @param w Z-velocity vector
   * @param i,j,k Grid indices
   * @param t Current time
   * @return Computed F function value
   */
  Real functionF_u(const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &w, int i, int j, int k, Real t);

  /**
   * @brief Computes the F function for v-velocity component
   */
  Real functionF_v(const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &w, int i, int j, int k, Real t);

  /**
   * @brief Computes the F function for w-velocity component
   */
  Real functionF_w(const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &w, int i, int j, int k, Real t);


  /**
   * @brief Manufactured solution function for u-velocity
   */
  Real functionG_u(int i, int j, int k, Real t);


  /**
   * @brief Manufactured solution function for v-velocity
   */
  Real functionG_v(int i, int j, int k, Real t);

  /**
   * @brief Manufactured solution function for w-velocity
   */
  Real functionG_w(int i, int j, int k, Real t);

  /**
   * @brief Main solution loop
   */
  void solve();

  /**
   * @brief Advances solution by one time step
   * @param time Current simulation time
   */
  void solve_time_step(Real time);

  /**
   * @brief Computes L2 error for X-velocity component
   * @param t Current time
   * @return L2 error value
   */
  Real error_comp_X(const Real t);

  /**
   * @brief Computes L2 error for Y-velocity component
   */
  Real error_comp_Y(const Real t);

  /**
   * @brief Computes L2 error for Z-velocity component
   */
  Real error_comp_Z(const Real t);

  /**
   * @brief Computes total L2 error
   */
  Real L2_error(const Real t);

  /**
   * @brief Parses input configuration file
   * @param input_file Path to input file
   */
  void parse_input(const std::string& input_file);

  // TODO write the output file.
  void output();

private:

  /// @brief MPI rank of current process
  int rank, size;

  /// @brief MPI Cartesian communicator
  MPI_Comm cart_comm;

  /// @brief Domain decomposition dimensions
  int dims[2];

  /// @brief Periodic boundary conditions flags
  int periods[2] = {1, 1};

  /// @brief Boundary conditions handler
  Boundary boundary;

  /// @brief Exact solution computer
  ExactSolution exact_solution;

  /// @brief Intermediate solution vectors
  std::vector<Real> Y2_x{}, Y2_y{}, Y2_z{};
  std::vector<Real> Y3_x{}, Y3_y{}, Y3_z{};

  /// @brief Input/output file paths
  std::string input_file;  // input file.
  std::string output_file; // output file.

  /// @brief Boundary flags for domain decomposition
  int lbx = 0, rbx = 0, lby = 0, rby = 0;

  /// @brief Brief to handle boundary cores()
  int firstX=0,firstY=0,lastX=0,lastY=0;

  /// @brief Process coordinates and neighbors in cartesian grid
  int coords[2];
  int neighbors[4];

  /// @brief Local grid data for each direction
  std::vector<Real> grid_loc_x{},grid_loc_y{}, grid_loc_z{};

  /// @brief MPI datatypes for face communication
  MPI_Datatype MPI_face_x_x, MPI_face_y_x;
  MPI_Datatype MPI_face_x_y, MPI_face_y_y;
  MPI_Datatype MPI_face_x_z, MPI_face_y_z;

  /// @brief MPI request handles for non-blocking communication
  MPI_Request reqs[4];

  /// @brief Grid dimensions and decomposition parameters
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
