/**
 * @file core.hpp
 * @brief Core implementation of the incompressible Navier-Stokes solver.
 *
 * This file contains the main solver class for incompressible Navier-Stokes equations.
 * It handles grid management, parallel computation, boundary conditions, and time-stepping
 * solutions using MPI for parallel processing.
 */

#ifndef CORE_HPP
#define CORE_HPP

#include "real.hpp"
#include "utils.hpp"
#include "boundary.hpp"
#include "grid.hpp"
#include "constants.hpp"
#include "C2Decomp.hpp"
#include <string>
#include <cmath>
#include <filesystem>
#include <mpi.h>
#include <fftw3.h>
#include "poissonSolver.hpp"
#include <math.h>
#include <fstream>
#include <memory>

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
   * @brief Constructor for the IcoNS solver.
   * 
   * @param cart_comm MPI Cartesian communicator
   * @param input_file Path to input configuration file
   * @param rank MPI process rank
   * @param size Total number of MPI processes
   */
  IcoNS(MPI_Comm cart_comm, const std::string &input_file, int rank, int size)
      : cart_comm(cart_comm),
        input_file(input_file),
        rank(rank),
        size(size)
  {
    parse_input(input_file);

    set2Decomp();

    setPoissonSolver();

    setBoundaryConditions();

    setParallelization();

    boundary.initializeBoundary(
        dim_x_x, dim_y_x, dim_x_y, dim_y_y,
        dim_x_z, dim_y_z, dim_z, dim_z_z,
        newDimX_x, newDimY_x, newDimX_y, newDimY_y,
        newDimX_z, newDimY_z,
        c2d->xSize);

  }

  /**
   * @brief Parses input configuration file
   * 
   * @param input_file Path to input file
   */
  void parse_input(const std::string& input_file);

  /**
   * @brief Sets up 2decomp library for domain decomposition.
   */
  void set2Decomp();

  /**
   * @brief Sets up Poisson solver for pressure calculation 
   *        depending on the test case.
   */
  void setPoissonSolver();

  /**
   * @brief Sets up boundary conditions for the simulation.
   */
  void setBoundaryConditions();

  /**
   * @brief Sets up parallel communication patterns and domain decomposition.
   */
  void setParallelization();

  /**
   * @brief Exchanges ghost cell data between neighboring processes
   * 
   * @param grid_loc Local grid data
   * @param newDimX New X dimension after decomposition
   * @param newDimY New Y dimension after decomposition
   * @param dim_z Z dimension
   * @param MPI_face_x MPI datatype for X-direction face exchange
   * @param MPI_face_y MPI datatype for Y-direction face exchange
   * @param sameX Only for periodic boundary, skip identical datas during communication(set to 0 only for communication x grids)
   * @param sameY Only for periodic boundary, skip identical datas during communication(set to 0 only for communication y grids)
   */
  void exchangeData(std::vector<Real> &grid_loc, int newDimX, int newDimY, int dim_z, 
                    MPI_Datatype MPI_face_x, MPI_Datatype MPI_face_y,int sameX, int sameY);

  /**
   * @brief Pressure copy for halo exchange.
   * 
   * @param p Pressure vector
   * @param halo pressure vector containing halo points
   */
  void copyPressureToHalo(Real* p, std::vector<Real> &halo);

  /**
   * @brief Computes the F function for u-velocity component.
   * 
   * @param u X-velocity vector
   * @param v Y-velocity vector
   * @param w Z-velocity vector
   * @param i,j,k Grid indices
   * @param t Current time
   * @return Computed F function value
   */
  inline Real functionF_u(const std::vector<Real> &u, const std::vector<Real> &v, 
                          const std::vector<Real> &w, int i, int j, int k, Real t);

  /**
   * @brief Computes the F function for v-velocity component.
   * 
   * @see functionF_u for parameter and return details.
   */
  inline Real functionF_v(const std::vector<Real> &u, const std::vector<Real> &v, 
                          const std::vector<Real> &w, int i, int j, int k, Real t);

  /**
   * @brief Computes the F function for w-velocity component.
   * 
   * @see functionF_u for parameter and return details.
   */
  inline Real functionF_w(const std::vector<Real> &u, const std::vector<Real> &v, 
                          const std::vector<Real> &w, int i, int j, int k, Real t);

  /**
   * @brief Manufactured solution function for u-velocity.
   */
  Real functionG_u(int i, int j, int k, Real t);

  /**
   * @brief Manufactured solution function for v-velocity.
   */
  Real functionG_v(int i, int j, int k, Real t);

  /**
   * @brief Manufactured solution function for w-velocity.
   */
  Real functionG_w(int i, int j, int k, Real t);

  /**
   * @brief Main solution loop.
   */
  void solve();

  /**
   * @brief Advances solution by one time step.
   * 
   * @param time Current simulation time
   */
  void solve_time_step(Real time);

  /**
   * @brief Computes L2 error for X-velocity component.
   * 
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
   * @brief Computes L2 error for pressure component
   */
  Real error_comp_P(const Real t);

  /**
   * @brief Computes total L2 error
   */
  void L2_error(const Real t);

  /**
   * @brief Output functions at final timestep
   */
  void output();
  void output_x();
  void output_y();
  void output_z();
  void output_profile();
  void cleaner();
  MPI_Status status;

private:

  // Test case number.
  int testCase;

  // MPI Cartesian communicator
  MPI_Comm cart_comm;

  // Domain decomposition dimensions
  int dims[2];

  // Periodic boundary conditions flags
  int periods[2] = {1, 1};

  // Boundary conditions handler
  Boundary boundary;

  // Exact solution object
  ExactSolution exact_solution;

  // Grid data structures
  Grid grid;
  Y2Grid y2Grid;
  Y3Grid y3Grid;

  // 2decomp library object
  C2Decomp *c2d;

  // Poisson solver object
  PoissonSolver *poissonSolver;

  // Arrays to store the dimensions of the grid in each direction
  int xSize[3], ySize[3], zSize[3];

  // Intermediate solution vector.
  std::vector<Real> halo_phi{};
  Real* Y2_p;

  // Periodicity for 2decomp
  bool periodss[3] = {false, true, true};

  // Input file.
  const std::string input_file;  

  // MPI rank of current process.
  int rank, size;

  // Boundary flags for domain decomposition.
  int lbx = 0, rbx = 0, lby = 0, rby = 0,lbz=0,rbz=0;

  // Brief to handle boundary cores.
  int firstX=0,firstY=0,lastX=0,lastY=0;

  // Process coordinates and neighbors in cartesian grid.
  int coords[2];
  int neighbors[4];

  // MPI datatypes for face communication.
  MPI_Datatype MPI_face_x_x, MPI_face_y_x;
  MPI_Datatype MPI_face_x_y, MPI_face_y_y;
  MPI_Datatype MPI_face_x_z, MPI_face_y_z;
  MPI_Datatype MPI_face_x_p, MPI_face_y_p;

  // MPI request handles for non-blocking communication.
  MPI_Request reqs[4];

  // Grid dimensions and decomposition parameters.
  int dim_x_x, dim_y_x;
  int dim_x_y, dim_y_y, dim_y_z;
  int dim_z, dim_x_z, dim_z_z;

  // New dimension parameters for each mesh after decomposition.
  int newDimX_x, newDimY_x;
  int newDimX_y, newDimY_y;
  int newDimX_z, newDimY_z;

  // Offsets for each mesh.
  int offset_x_x,offset_y_x, offset_x_y, offset_y_y, offset_x_z, offset_y_z;
  int resx = 0,resy = 0;

  // functions to access easier to each grid.
  inline int getx(int i, int j, int k) 
  {
    if(k==-1)
    {
      k=dim_z-2;
    }
    if(k==dim_z)
    {
      k=1;
    }
     return i * newDimY_x * dim_z + j * dim_z + k; 
  }

  inline int gety(int i, int j, int k) 
  {
    if (k == -1)
    {
      k = dim_z - 2;
    }
    if (k == dim_z)
    {
      k = 1;
    }
     return i * newDimY_y * dim_z + j * dim_z + k; 
  }

  inline int getz(int i, int j, int k) { return i * newDimY_z * dim_z_z + j * dim_z_z + (k+(dim_z_z))%(dim_z_z); }

  inline int getp(int i, int j, int k) 
  {
    if (k == -1)
    {
      k = xSize[0] - 2;
    }
    if (k == xSize[0])
    {
      k = 1;
    }
     return i * xSize[1] * xSize[0] + j * xSize[0] + k; 
  }

  inline int getHaloP(int i, int j, int k) 
  {
    if (k == -1)
    {
      k = xSize[0] - 2;
    }
    if (k == xSize[0])
    {
      k = 1;
    }
     return i * (xSize[1] + 2) * xSize[0] + j * xSize[0] + k; 
  }
};

#endif // CORE_HPP

