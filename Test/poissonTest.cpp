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
    // std::array<Real, (NX+1) * (NY+1) * (NZ+1)> F;
    // std::array<Real, (NX+1) * (NY+1) * (NZ+1)> sol;
    
    for (size_t i=0; i < NX; i++){
        for (size_t j=0; j < NY; j++){
            for (size_t k=0; k < NZ; k++){
                F[i * (NY) * (NZ) + j * (NZ) + k] = ((DX/LX) * (DX/LX) + (DY/LY) * (DY/LY) + (DZ/LZ) * (DZ/LZ)) *
                                                          4 * M_PI * M_PI * 
                                                          std::cos(i * DX * 2 * M_PI/LX) *
                                                          std::cos(j * DY * 2 * M_PI/LY) * 
                                                          std::cos(k * DZ * 2 * M_PI/LZ);

                sol[i * (NY) * (NZ) + j * (NZ) + k] = -std::cos(i * DX * 2 * M_PI/LX) *
                                                            std::cos(j * DY * 2 * M_PI/LY) *
                                                            std::cos(k * DZ * 2 * M_PI/LZ);
            }   
        }
    }

    fftw_complex* out = fftw_alloc_complex(NX * NY * (NZ/2 + 1));

    poissonSolver.solveDirichletPoisson(F, out);

    fftw_free(out);

    double error = 0;

    for (size_t i=0; i < NX; i++){
        for (size_t j=0; j < NY; j++){
            for (size_t k=0; k < NZ; k++){
                error = error + ((sol[i * (NY) * (NZ) + j * (NZ) + k] - F[i * (NY) * (NZ) + j * (NZ) + k]) *
                                 (sol[i * (NY) * (NZ) + j * (NZ) + k] - F[i * (NY) * (NZ) + j * (NZ) + k]) *
                                 (DX * DY * DZ));
            }   
        }
    }

    std::cout << sqrt(error) << std::endl;

    // for (size_t i=0; i < NX+1; i++){
    //     for (size_t j=0; j < NY+1; j++){
    //         for (size_t k=0; k < NZ+1; k++){
    //             F[i * (NY+1) * (NZ+1) + j * (NZ+1) + k] = ((DX/LX) * (DX/LX) + (DY/LY) * (DY/LY) + (DZ/LZ) * (DZ/LZ)) *
    //                                                       4 * M_PI * M_PI * 
    //                                                       std::cos(i * DX * 2 * M_PI/LX) *
    //                                                       std::cos(j * DY * 2 * M_PI/LY) * 
    //                                                       std::cos(k * DZ * 2 * M_PI/LZ);

    //             sol[i * (NY+1) * (NZ+1) + j * (NZ+1) + k] = -std::cos(i * DX * 2 * M_PI/LX) *
    //                                                         std::cos(j * DY * 2 * M_PI/LY) *
    //                                                         std::cos(k * DZ * 2 * M_PI/LZ);
    //         }   
    //     }
    // }

    // poissonSolver.solveNeumannPoisson(F);

    // double error = 0;

    // for (size_t i=0; i < NX+1; i++){
    //     for (size_t j=0; j < NY+1; j++){
    //         for (size_t k=0; k < NZ+1; k++){
    //             error = error + ((sol[i * (NY+1) * (NZ+1) + j * (NZ+1) + k] - F[i * (NY+1) * (NZ+1) + j * (NZ+1) + k]) *
    //                              (sol[i * (NY+1) * (NZ+1) + j * (NZ+1) + k] - F[i * (NY+1) * (NZ+1) + j * (NZ+1) + k]) *
    //                              (DX * DY * DZ));
    //         }   
    //     }
    // }

    // std::cout << sqrt(error) << std::endl;

    return 0;
}
