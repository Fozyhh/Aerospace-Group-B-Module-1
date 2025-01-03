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

    C2Decomp *c2d;
    int xSize[3], ySize[3], zSize[3];

public: 
    PoissonSolver(const bool periodicX, const bool periodicY, const bool periodicZ, 
                  C2Decomp *c2d)
    : periodicX(periodicX),
      periodicY(periodicY),
      periodicZ(periodicZ),
      c2d(c2d),
      xSize{c2d->xSize[0], c2d->xSize[1], c2d->xSize[2]},
      ySize{c2d->ySize[0], c2d->ySize[1], c2d->ySize[2]},
      zSize{c2d->zSize[0], c2d->zSize[1], c2d->zSize[2]}
    {}

    //TODO: array brutti, 2decomp per la separazione
    void solveDirichletPoisson(std::vector<Real>& F_dP, fftw_complex *FD);
    void solveNeumannPoisson(double* F);
};

#endif

// We have to decide how to pass different types of boundaries (Periodic or Neumann). For now we will use a boolean