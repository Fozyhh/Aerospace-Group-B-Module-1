/**
 * @file core.hpp
 * @brief Core solver for incompressible Navier-Stokes equations
 * @details Contains the main solver class for handling the numerical solution
 *          of incompressible Navier-Stokes equations including the grid, the time step,
 *          the final time, the initial condition, the boundary conditions, the source term,
 *          the exact solution, the numerical solution, the error, and the output file.
 */

#ifndef CORE_HPP
#define CORE_HPP

#include "utils.hpp"
#include "boundary.hpp"
#include "grid.hpp"
#include <string>
#include <cmath>

/**
 * @class IcoNS
 * @brief Incompressible Navier-Stokes solver
 * @details Handles the numerical solution of incompressible Navier-Stokes equations
 *          including grid, time integration, boundary conditions, ...
 */
class IcoNS
{
public:
  
  /**
   * @brief Constructor for the IcoNS solver
   * @param lx Domain length in x direction
   * @param ly Domain length in y direction
   * @param lz Domain length in z direction
   * @param nx Number of cells in x direction
   * @param ny Number of cells in y direction
   * @param nz Number of cells in z direction
   * @param dt Time step size
   * @param T Final simulation time
   * @param Re Reynolds number
   * @param input_file Path to input file
   * @param output_file Path to output file
   */
  IcoNS(const double lx, const double ly, const double lz,
        const unsigned int nx, const unsigned int ny, const unsigned int nz,
        const double dt, const double T, const double Re,
        const std::string &input_file, const std::string &output_file)
      : grid(nx, ny, nz),
        dt(dt),
        T(T),
        Re(Re),
        lx(lx),
        ly(ly),
        lz(lz),
        nx(nx),
        ny(ny),
        nz(nz),
        dx(lx / nx),
        dy(ly / ny),
        dz(lz / nz),
        input_file(input_file),
        output_file(output_file),
        boundary(&grid,dx,dy,dz)
  {}

  /**
   * @brief Initialize the grid and simulation parameters
   */
  void preprocessing(/*std::string &input_file*/);

  /**
   * @brief Compute the source term F
   * @param u X-velocity component
   * @param v Y-velocity component
   * @param w Z-velocity component
   * @param i X-index
   * @param j Y-index
   * @param k Z-index
   * @param t Current time
   * @return Vector containing source term components
   */
  std::vector<double> functionF(const std::vector<double> &u, const std::vector<double> &v, const std::vector<double> &w, 
                                size_t i, size_t j, size_t k, double t);
 
  /**
   * @brief Compute the source term G
   * @param i X-index
   * @param j Y-index
   * @param k Z-index
   * @param t Current time
   * @return Vector containing source term components
   */
  std::vector<double> functionG(size_t i, size_t j, size_t k, double t);

  /**
   * @brief Apply boundary conditions at given time
   * @param time Current simulation time
   */
  void apply_boundary_conditions(double time);
  
  /**
   * @brief Solve a single time step
   * @param time Current simulation time
   */
  void solve_time_step( double time );
  
  /**
   * @brief Solve the complete problem and save output
   */
  void solve();

  /**
   * @brief Alternative source term computation
   * @param u X-velocity component
   * @param v Y-velocity component
   * @param w Z-velocity component
   * @param i X-index
   * @param j Y-index
   * @param k Z-index
   * @return Vector containing source term components
   */
  std::vector<double> functionF(const std::vector<double> &u, const std::vector<double> &v, const std::vector<double> &w,
                                size_t i, size_t j, size_t k);
  
  /**
   * @brief Alternative time step solver
   */
  void solve_time_step();

  /**
   * @brief Compute error in X component
   * @param t Current time
   * @return Error value
   */
  double error_comp_X(const double t);

  /**
   * @brief Compute error in Y component
   * @param t Current time
   * @return Error value
   */
  double error_comp_Y(const double t);
  
  /**
   * @brief Compute error in Z component
   * @param t Current time
   * @return Error value
   */ 
  double error_comp_Z(const double t);
  
  /**
   * @brief Compute the L2 norm of the error
   * @param t Current time
   * @return L2 norm value
   */
  double L2_error(const double t);

  /**
   * @brief Write results to output file
   */
  void output(); // write the output file.


private:
    Grid grid;                     ///< grid of the domain.
    const double dt;               ///< time step.
    const double T;                ///< final time.
    const double Re;               ///< Reynolds number.
    const unsigned int lx, ly, lz; ///< lengths of edges of the domain.
    const unsigned int nx, ny, nz; ///< number of cells in the x,y,z directions.
    const double dx, dy, dz;       ///< cell sizes in the x,y,z directions.
    Boundary boundary;
    ExactSolution exact_solution;  ///< exact solution.
    std::string input_file;        ///< input file.
    std::string output_file;       ///< output file.

};

#endif // CORE_HPP
