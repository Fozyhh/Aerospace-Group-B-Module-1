// class representing the problem to be solved. It contains the grid, the time step, the final time, the initial condition, the boundary conditions, the source term, the exact solution, the numerical solution, the error, and the output file.
// the boundary conditions, the source term, the exact solution, the numerical solver.

#ifndef CORE_HPP
#define CORE_HPP

#include "grid.hpp"

class IcoNS
{
public:
    class ExactSolution
    {
    };

    IcoNS(const unsigned int nx, const unsigned int ny, const unsigned int nz, const double dt, const double T, const double Re, const std::string &input_file, const std::string &output_file)
        : grid(nx, ny, nz),
          dt(dt),
          T(T),
          Re(Re),
          nx(nx),
          ny(ny),
          nz(nz),
          input_file(input_file),
          output_file(output_file)
    {
    }

    void preprocessing(string &input_file); // grid initialization.

    void solve(); // solve the problem saving the ouput.

    void solve_time_step(); // solve a time step.

    void output(); // write the output file.

private:
    Grid grid;                     // grid of the domain.
    const double dt;               // time step.
    const double T;                // final time.
    const double Re;               // Reynolds number.
    const unsigned int nx, ny, nz; // number of cells in the x,y,z directions.
    ExactSolution exact_solution;  // exact solution.
    std::string input_file;        // input file.
    std::string output_file;       // output file.
};
