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
  IcoNS(const Real lx, const Real ly, const Real lz,
        const unsigned int nx, const unsigned int ny, const unsigned int nz,
        const double dt, const double T, const Real Re,
        const std::string &input_file, const std::string &output_file)
      : grid(nx, ny, nz),
        boundary(nx, ny, nz, lx/nx, ly/ny, lz/nz),
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
        exact_solution(lx/nx, ly/ny, lz/nz)
  {
  }

  void preprocessing(/*std::string &input_file*/); // grid initialization.

  Real functionF_u(const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &w, size_t i, size_t j, size_t k, Real t); // compute the source term.
  Real functionF_v(const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &w, size_t i, size_t j, size_t k, Real t); // compute the source term.
  Real functionF_w(const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &w, size_t i, size_t j, size_t k, Real t); // compute the source term.
  Real functionG_u(size_t i, size_t j, size_t k, Real t);                                                                                           // compute the source term.
  Real functionG_v(size_t i, size_t j, size_t k, Real t);                                                                                           // compute the source term.
  Real functionG_w(size_t i, size_t j, size_t k, Real t);                                                                                           // compute the source term.

  void apply_boundary_conditions(double time); // apply the boundary conditions.
  void solve_time_step(double time);           // solve a time step.
  void solve();                                // solve the problem saving the ouput.

  Real error_comp_X(const Real t);
  Real error_comp_Y(const Real t);
  Real error_comp_Z(const Real t);
  Real L2_error(const Real t); // compute the L2 norm

  void output(); // write the output file.

private:
  Grid grid; // grid of the domain.
  Boundary boundary;
  ExactSolution exact_solution;  // exact solution.
  const double dt;               // time step.
  const double T;                // final time.
  const Real Re;               // Reynolds number.
  const Real lx, ly, lz; // lengths of edges of the domain.
  const unsigned int nx, ny, nz; // number of cells in the x,y,z directions.
  const Real dx, dy, dz;       // cell sizes in the x,y,z directions.
  std::string input_file;        // input file.
  std::string output_file;       // output file.
};

#endif // CORE_HPP
