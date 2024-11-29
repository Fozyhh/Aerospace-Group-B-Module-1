#include <iostream>
#include <math.h>
#include "../includes/poissonSolver.hpp"
#include <fftw3.h>


int main()
{
    const bool a = true;
    const bool b = true;
    const bool c = true;

    PoissonSolver poissonSolver(a,b,c);

    std::array<Real, NX * NY * NZ> F;
    std::array<Real, NX * NY * NZ> sol;

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

    // F and sol are OK

    poissonSolver.solvePoisson(F, out);

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

    return 0;
}
