#include <iostream>
#include <math.h>
#include "poissonSolver.hpp"
#include <fftw3.h>


int main()
{
    const bool a = true;
    const bool b = true;
    const bool c = true;

    PoissonSolver poissonSolver(a,b,c);

    std::array<Real, NX * NY * NZ> F;
    std::array<Real, NX * NY * NZ> sol;
    
    // for (size_t i=0; i < NX; i++){
    //     for (size_t j=0; j < NY; j++){
    //         for (size_t k=0; k < NZ; k++){
    //             F[i * (NY) * (NZ) + j * (NZ) + k] = ((DX/LX) * (DX/LX) + (DY/LY) * (DY/LY) + (DZ/LZ) * (DZ/LZ)) *
    //                                                       4 * M_PI * M_PI * 
    //                                                       std::cos(i * DX * 2 * M_PI/LX) *
    //                                                       std::cos(j * DY * 2 * M_PI/LY) * 
    //                                                       std::cos(k * DZ * 2 * M_PI/LZ);

    //             sol[i * (NY) * (NZ) + j * (NZ) + k] = -std::cos(i * DX * 2 * M_PI/LX) *
    //                                                         std::cos(j * DY * 2 * M_PI/LY) *
    //                                                         std::cos(k * DZ * 2 * M_PI/LZ);
    //         }   
    //     }
    // }

    // fftw_complex* out = fftw_alloc_complex(NX * NY * (NZ/2 + 1));

    // // F and sol are OK

    // poissonSolver.solveDirichletPoisson(F, out);

    // fftw_free(out);

    // double error = 0;

    // for (size_t i=0; i < NX; i++){
    //     for (size_t j=0; j < NY; j++){
    //         for (size_t k=0; k < NZ; k++){
    //             error = error + ((sol[i * (NY) * (NZ) + j * (NZ) + k] - F[i * (NY) * (NZ) + j * (NZ) + k]) *
    //                              (sol[i * (NY) * (NZ) + j * (NZ) + k] - F[i * (NY) * (NZ) + j * (NZ) + k]) *
    //                              (DX * DY * DZ));
    //         }   
    //     }
    // }

    // std::cout << sqrt(error) << std::endl;

    // using DX, DY, DZ gives first order accuracy
    double dx = LX/(NX-1);
    double dy = LY/(NY-1);
    double dz = LZ/(NZ-1);
    // using dx, dy, dz gives second order accuracy (starnge stuff, is it fine?)

    for (size_t i=0; i < NX; i++){
        for (size_t j=0; j < NY; j++){
            for (size_t k=0; k < NZ; k++){
                F[i * (NY) * (NZ) + j * (NZ) + k] = ((dx/LX) * (dx/LX) + (dy/LY) * (dy/LY) + (dz/LZ) * (dz/LZ)) *
                                                          4 * M_PI * M_PI * 
                                                          std::cos(i * dx * 2 * M_PI/LX) *
                                                          std::cos(j * dy * 2 * M_PI/LY) * 
                                                          std::cos(k * dz * 2 * M_PI/LZ);

                sol[i * (NY) * (NZ) + j * (NZ) + k] = -std::cos(i * dx * 2 * M_PI/LX) *
                                                            std::cos(j * dy * 2 * M_PI/LY) *
                                                            std::cos(k * dz * 2 * M_PI/LZ);
            }   
        }
    }

    poissonSolver.solveNeumannPoisson(F);

    double error = 0;

    for (size_t i=0; i < NX; i++){
        for (size_t j=0; j < NY; j++){
            for (size_t k=0; k < NZ; k++){
                error = error + ((sol[i * (NY) * (NZ) + j * (NZ) + k] - F[i * (NY) * (NZ) + j * (NZ) + k]) *
                                 (sol[i * (NY) * (NZ) + j * (NZ) + k] - F[i * (NY) * (NZ) + j * (NZ) + k]) *
                                 (dx * dy * dz));
            }   
        }
    }

    std::cout << sqrt(error) << std::endl;

    return 0;
}
