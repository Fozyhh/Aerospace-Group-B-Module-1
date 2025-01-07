#ifndef POISSONSOLVER_HPP
#define POISSONSOLVER_HPP

#include "utils.hpp"
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

private:
    /// @brief Flag for periodic boundary condition along the x-axis
    const bool periodicX;
    /// @brief Flag for periodic boundary condition along the y-axis
    const bool periodicY;
    /// @brief Flag for periodic boundary condition along the z-axis
    const bool periodicZ;
    
    /// @brief Pointer to the C2Decomp object for parallel decomposition
    C2Decomp *c2d;

    /// @brief Arrays holding grid size information for each view
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

      /**
     * @brief Solves Poisson's equation with Dirichlet boundary conditions.
     * 
     * @param F_dP A vector containing the source term.
     * @param FD A pointer to an FFTW complex array for frequency-domain calculations.
     */
    void solveDirichletPoisson(std::vector<Real>& F_dP, fftw_complex *FD);


    /**
     * @brief Solves Poisson's equation with Neumann boundary conditions.
     * 
     * @param F A pointer to the grid data for Neumann boundary conditions.
     */
    void solveNeumannPoisson(double* F);
};

#endif