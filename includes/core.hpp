/// Class representing the problem to be solved. It contains the grid, the time step, the final time, the initial condition, the boundary conditions, the source term, the exact solution, the numerical solution, the error, and the output file.
/// the boundary conditions, the source term, the exact solution, the numerical solver.

#ifndef CORE_HPP
#define CORE_HPP

#include "utils.hpp"
#include "boundary.hpp"
#include "grid.hpp"
#include <fftw3.h>
#include <string>
#include <cmath>
#include <filesystem>
#define PERIODIC
//#define DIRICHELET

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
  Real functionG_u(size_t i, size_t j, size_t k, Real t);                                                                                                                                                              // compute the source term.
  Real functionG_v(size_t i, size_t j, size_t k, Real t);                                                                                                                                                              // compute the source term.
  Real functionG_w(size_t i, size_t j, size_t k, Real t);                                                                                                                                                              // compute the source term.

  void solve();                    // solve the problem saving the ouput.
  void solve_time_step(Real time); // solve a time step.

  Real error_comp_X(const Real t);
  Real error_comp_Y(const Real t);
  Real error_comp_Z(const Real t);
  Real error_comp_P(const Real t);
  Real l2_norm_x(const Real t);
  Real l2_norm_y(const Real t);
  Real l2_norm_z(const Real t);
  Real L2_error(const Real t); // compute the L2 norm

  void output_u(const std::string &filename, Grid &print); // write the output file.
  void output_v(const std::string &filename, Grid &print); // write the output file.
  void output_w(const std::string &filename, Grid &print); // write the output file.
  std::array<Real, (NX * (NY + 1) * (NZ + 1))> u_function(){return grid.u;};
  std::array<Real, ((NX + 1) * NY * (NZ + 1))> v_function(){return grid.v;};
  std::array<Real, ((NX + 1) * (NY + 1) * NZ)> w_function(){return grid.w;};
  fftw_complex* helper = fftw_alloc_complex(NX * NY * (NZ/2 + 1));
private:
  Grid grid; // grid of the domain.
  Boundary boundary;
  ExactSolution exact_solution;
  std::array<Real, (NX * (NY + 1) * (NZ + 1))> Y2_x{};
  std::array<Real, ((NX + 1) * NY * (NZ + 1))> Y2_y{};
  std::array<Real, ((NX + 1) * (NY + 1) * NZ)> Y2_z{};
  std::array<Real, ((NX) * (NY) * (NZ))> Y2_p{};
  //std::array<Real, ((NX) * (NY) * (NZ))> Sol_p{};
  std::array<Real, ((NX) * (NY) * (NZ))> Phi_p{};
  std::array<Real, (NX * (NY + 1) * (NZ + 1))> Y3_x{};
  std::array<Real, ((NX + 1) * NY * (NZ + 1))> Y3_y{};
  std::array<Real, ((NX + 1) * (NY + 1) * NZ)> Y3_z{};
  std::array<Real, ((NX) * (NY) * (NZ))> Y3_p{};
  std::string input_file;  // input file.
  std::string output_file; // output file.


#ifdef PERIODIC
    inline size_t indexingPeriodicx(size_t i, size_t j, size_t k) { return ((i + NX) % NX) * (NY + 1) * (NZ + 1) + ((j + NY) % NY) * (NZ + 1) + ((k + NZ) % NZ); };
    inline size_t indexingPeriodicy(size_t i, size_t j, size_t k) { return ((i + NX) % NX) * NY * (NZ + 1) + ((j + NY) % NY) * (NZ + 1) + ((k + NZ) % NZ); };
    inline size_t indexingPeriodicz(size_t i, size_t j, size_t k) { return ((i + NX) % NX) * (NY + 1) * NZ + ((j + NY) % NY) * NZ + ((k + NZ) % NZ); };
    inline size_t indexingPeriodicp(size_t i, size_t j, size_t k) { return ((i + NX) % NX) * (NY + 1) * (NZ + 1) + ((j + NY) % NY) * (NZ + 1) + ((k + NZ) % NZ); };
#endif
  
#ifdef DIRICHELET
    inline size_t indexingDiricheletx(size_t i, size_t j, size_t k) { return i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k; }
    inline size_t indexingDirichelety(size_t i, size_t j, size_t k) { return i * NY * (NZ + 1) + j * (NZ + 1) + k; }
    inline size_t indexingDiricheletz(size_t i, size_t j, size_t k) { return i * (NY + 1) * NZ + j * NZ + k; }
#endif
};

#endif // CORE_HPP
