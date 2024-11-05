/// Class representing the problem to be solved. It contains the grid, the time step, the final time, the initial condition, the boundary conditions, the source term, the exact solution, the numerical solution, the error, and the output file.
/// the boundary conditions, the source term, the exact solution, the numerical solver.

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
  IcoNS(const std::string &input_file, const std::string &output_file)
      : input_file(input_file),
        output_file(output_file)
  {
  }

  void preprocessing(/*std::string &input_file*/); // grid initialization.

  double functionF_u(const std::array<double, NX *(NY + 1) * (NZ + 1)> &u, const std::array<double, (NX + 1) * NY *(NZ + 1)> &v, const std::array<double, (NX + 1) * (NY + 1) * NZ> &w, size_t i, size_t j, size_t k, double t); // compute the source term.
  double functionF_v(const std::array<double, NX *(NY + 1) * (NZ + 1)> &u, const std::array<double, (NX + 1) * NY *(NZ + 1)> &v, const std::array<double, (NX + 1) * (NY + 1) * NZ> &w, size_t i, size_t j, size_t k, double t); // compute the source term.
  double functionF_w(const std::array<double, NX *(NY + 1) * (NZ + 1)> &u, const std::array<double, (NX + 1) * NY *(NZ + 1)> &v, const std::array<double, (NX + 1) * (NY + 1) * NZ> &w, size_t i, size_t j, size_t k, double t); // compute the source term.
  double functionG_u(size_t i, size_t j, size_t k, double t);                                                                                                                                                                    // compute the source term.
  double functionG_v(size_t i, size_t j, size_t k, double t);                                                                                                                                                                    // compute the source term.
  double functionG_w(size_t i, size_t j, size_t k, double t);                                                                                                                                                                    // compute the source term.

  void solve();                      // solve the problem saving the ouput.
  void solve_time_step(double time); // solve a time step.

  double error_comp_X(const double t);
  double error_comp_Y(const double t);
  double error_comp_Z(const double t);
  double L2_error(const double t); // compute the L2 norm

  void output(); // write the output file.

private:
  Grid grid; // grid of the domain.
  Boundary boundary;
  ExactSolution exact_solution;
  std::array<double, (NX * (NY + 1) * (NZ + 1))> Y2_x{};
  std::array<double, ((NX + 1) * NY * (NZ + 1))> Y2_y{};
  std::array<double, ((NX + 1) * (NY + 1) * NZ)> Y2_z{};
  std::array<double, (NX * (NY + 1) * (NZ + 1))> Y3_x{};
  std::array<double, ((NX + 1) * NY * (NZ + 1))> Y3_y{};
  std::array<double, ((NX + 1) * (NY + 1) * NZ)> Y3_z{};
  std::string input_file;  // input file.
  std::string output_file; // output file.
};

#endif // CORE_HPP
