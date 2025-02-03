#ifndef POISSONSOLVER_HPP
#define POISSONSOLVER_HPP

#include "utils.hpp"
#include "grid.hpp"
#include <cmath>
#include <numbers>
#include <complex>
#include <fftw3.h>
#include "../dependencies/2Decomp_C/C2Decomp.hpp"


/**
 * @class PoissonSolver
 * @brief Class handling solution of Poisson's equation 
 *
 * This class leverages the FFTW library for efficient Fourier transforms and integrates with the 
 * C2Decomp library for parallel decomposition. It provides methods for solving Poisson's equation.
 */

class PoissonSolver
{

protected:
    
    /// @brief Pointer to the C2Decomp object for parallel decomposition.
    C2Decomp *c2d;

    /// @brief Arrays holding grid size information for each view.
    int xSize[3], ySize[3], zSize[3];

public: 

    Real *pz;                // z-pencil
    Real *py;                // y_pencil
    FFTW_TYPE(plan) neumann;

    PoissonSolver(C2Decomp *c2d)
    : c2d(c2d),
      xSize{c2d->xSize[0], c2d->xSize[1], c2d->xSize[2]},
      ySize{c2d->ySize[0], c2d->ySize[1], c2d->ySize[2]},
      zSize{c2d->zSize[0], c2d->zSize[1], c2d->zSize[2]}
    {
      c2d->allocY(py);
      c2d->allocZ(pz);
    }

    /**
     * @brief Solves Poisson's equation. 
     * 
     * @param F A pointer to the grid data for pressure.
     */
    virtual void solvePoisson(Real* F) = 0;
};

/**
 * @brief Poisson solver for Neumann boundary conditions. (test1)
 */
class NeumannPoissonSolver : public PoissonSolver
{
public:
    NeumannPoissonSolver(C2Decomp *c2d): PoissonSolver(c2d)
    {}

    void solvePoisson(Real* F) override;
};

/**
 * @brief Poisson solver for Dirichlet boundary conditions. (test2)
 */
class DirichletPoissonSolver : public PoissonSolver
{
public:
    DirichletPoissonSolver(C2Decomp *c2d): PoissonSolver(c2d)
    {}

    void solvePoisson(Real* F) override;
};

#endif