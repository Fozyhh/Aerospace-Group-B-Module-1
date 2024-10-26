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
        const double dt, const double T, const double Re,
        const std::string &input_file, const std::string &output_file)
      : dt(dt),
        T(T),
        Re(Re),
        lx(lx),
        ly(ly),
        lz(lz),
        dx(lx / NX),
        dy(ly / NY),
        dz(lz / NZ),
        input_file(input_file),
        output_file(output_file)
  // boundary(&grid, dx, dy, dz)
  {
  }

  void preprocessing(/*std::string &input_file*/); // grid initialization.

  std::array<double, 3> functionF(const std::vector<std::vector<std::vector<double>>> &u,
                                  const std::vector<std::vector<std::vector<double>>> &v,
                                  const std::vector<std::vector<std::vector<double>>> &w,
                                  size_t i, size_t j, size_t k, double t); // compute the source term.

  std::array<double, 3> functionF(const std::array<std::array<std::array<double, NZ + 1>, NY + 1>, NX> &u,
                                  const std::array<std::array<std::array<double, NZ + 1>, NY>, NX + 1> &v,
                                  const std::array<std::array<std::array<double, NZ>, NY + 1>, NX + 1> &w,
                                  size_t i, size_t j, size_t k, double t);

  std::vector<double> functionG(size_t i, size_t j, size_t k, double t); // compute the source term.

  // void apply_boundary_conditions(double time); // apply the boundary conditions.
  void solve_time_step(double time); // solve a time step.
  void solve();                      // solve the problem saving the ouput.

  double error_comp_X(const double t);
  double error_comp_Y(const double t);
  double error_comp_Z(const double t);
  double L2_error(const double t); // compute the L2 norm

  void output(); // write the output file.

private:
  Grid grid;                     // grid of the domain.
  const double dt;               // time step.
  const double T;                // final time.
  const double Re;               // Reynolds number.
  const unsigned int lx, ly, lz; // lengths of edges of the domain.
  const double dx, dy, dz;       // cell sizes in the x,y,z directions.
  // Boundary boundary;
  ExactSolution exact_solution; // exact solution.
  std::string input_file;       // input file.
  std::string output_file;      // output file.
};

#endif // CORE_HPP
