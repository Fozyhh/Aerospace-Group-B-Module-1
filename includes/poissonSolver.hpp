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
    std::array<Real, NX * NY * NZ> eigenvalues;


public: 
    PoissonSolver(const bool periodicX, const bool periodicY, const bool periodicZ)
    : periodicX(periodicX),
      periodicY(periodicY),
      periodicZ(periodicZ)
    {
        // Here we will build the eigenvalues array -> Each eigenvalue is constructed considering the corresponding boolean of the direction
    }

    void solvePoisson(std::array<Real, NX * NY * NZ>& F, fftw_complex *FD);
    void solveNeumannPoisson(std::array<Real, (NX+1) * (NY+1) * (NZ+1)>& F);
};

#endif

// We have to decide how to pass different types of boundaries (Periodic or Neumann). For now we will use a boolean