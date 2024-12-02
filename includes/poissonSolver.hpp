#ifndef POISSONSOLVER_HPP
#define POISSONSOLVER_HPP

#include "utils.hpp"
#include <cmath>
#include <numbers>
#include <complex>
#include <fftw3.h>
#include "C2Decomp.hpp"

class PoissonSolver
{

private:
    const bool periodicX;
    const bool periodicY;
    const bool periodicZ;

public: 
    PoissonSolver(const bool periodicX, const bool periodicY, const bool periodicZ)
    : periodicX(periodicX),
      periodicY(periodicY),
      periodicZ(periodicZ)
    {}

    void solveDirichletPoisson(std::array<Real, (NX) * (NY) * (NZ)>& F_dP, fftw_complex *FD);
    void solveNeumannPoisson(std::array<Real, (NX+1) * (NY+1) * (NZ+1)>& F);
};

#endif

// We have to decide how to pass different types of boundaries (Periodic or Neumann). For now we will use a boolean