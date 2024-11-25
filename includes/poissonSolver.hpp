#ifndef POISSONSOLVER_HPP
#define POISSONSOLVER_HPP

#include "utils.hpp"

class PoissonSolver
{

private:
    const bool periodicX;
    const bool periodicY;
    const bool periodicZ;
    std::array<Real, (NX+1)*(NY+1)*(NZ+1)> eigenvalues;


public: 
    PoissonSolver(const bool periodicX, const bool periodicY, const bool periodicZ)
    : periodicX(periodicX),
      periodicY(periodicY),
      periodicZ(periodicZ)
    {
        // Here we will build the eigenvalues array -> Each eigenvalue is constructed considering the corresponding boolean of the direction
    }

    void solvePoisson(std::array<Real, (NX+1)*(NY+1)*(NZ+1)>& F, std::array<Real, (NX+1)*(NY+1)*(NZ+1)>& dP);
};

#endif

// We have to decide how to pass different types of boundaries (Periodic or Neumann). For now we will use a boolean