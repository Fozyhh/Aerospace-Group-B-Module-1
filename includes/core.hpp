// class representing the problem to be solved. It contains the grid, the time step, the final time, the initial condition, the boundary conditions, the source term, the exact solution, the numerical solution, the error, and the output file.
// the boundary conditions, the source term, the exact solution, the numerical solver.

#ifndef CORE_HPP
#define CORE_HPP

#include "utils.hpp"
#include "boundary.hpp"
#include "grid.hpp"
#include <string>
#include <cmath>

class IcoNS
{
public:
  IcoNS(const double lx, const double ly, const double lz,
        const unsigned int nx, const unsigned int ny, const unsigned int nz,
        const double dt, const double T, const double Re,
        const std::string &input_file, const std::string &output_file)
      : grid(nx, ny, nz),
        boundary(grid, lx/nx, ly/ny, lz/nz),
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
        output_file(output_file)
  {
  }

  void preprocessing(/*std::string &input_file*/); // grid initialization.

  std::vector<double> functionF(const std::vector<double> &u, const std::vector<double> &v, const std::vector<double> &w, size_t i, size_t j, size_t k, double t); // compute the source term.
  std::vector<double> functionG(size_t i, size_t j, size_t k, double t);                                                                                           // compute the source term.

  void apply_boundary_conditions(double time); // apply the boundary conditions.
  void solve_time_step(double time);           // solve a time step.
  void solve();                                // solve the problem saving the ouput.

  std::vector<double> functionF(const std::vector<double> &u, const std::vector<double> &v, const std::vector<double> &w,
                                size_t i, size_t j, size_t k);

  void solve_time_step(); // solve a time step.

  double error_comp_X(const double t);
  double error_comp_Y(const double t);
  double error_comp_Z(const double t);
  double L2_error(const double t); // compute the L2 norm

  void output(); // write the output file.

private:
  Grid grid; // grid of the domain.
  Boundary boundary;
  ExactSolution exact_solution;  // exact solution.
  const double dt;               // time step.
  const double T;                // final time.
  const double Re;               // Reynolds number.
  const unsigned int lx, ly, lz; // lengths of edges of the domain.
  const unsigned int nx, ny, nz; // number of cells in the x,y,z directions.
  const double dx, dy, dz;       // cell sizes in the x,y,z directions.
  std::string input_file;        // input file.
  std::string output_file;       // output file.
};

#endif // CORE_HPP
