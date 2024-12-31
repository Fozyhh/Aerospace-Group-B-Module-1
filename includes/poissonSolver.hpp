#ifndef POISSONSOLVER_HPP
#define POISSONSOLVER_HPP

#include "utils.hpp"
#include <cmath>
#include <numbers>
#include <complex>
#include <fftw3.h>
#include "../dependencies/2Decomp_C/C2Decomp.hpp"

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

    //TODO: array brutti, 2decomp per la separazione
    void solveDirichletPoisson(std::vector<Real>& F_dP, fftw_complex *FD);
    void solveNeumannPoisson(std::vector<Real>& F);
};

#endif

// We have to decide how to pass different types of boundaries (Periodic or Neumann). For now we will use a boolean