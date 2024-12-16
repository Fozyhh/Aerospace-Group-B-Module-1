#include <iostream>
#include <math.h>
#include "poissonSolver.hpp"
#include <fftw3.h>
#include <mpi.h>


int main(int argc, char *argv[])
{
    int ierr, totRank, mpiRank;

    //Initialize MPI
    ierr = MPI_Init( &argc, &argv);

    //Get the number of processes
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &totRank);

    //Get the local rank
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    const bool a = true;
    const bool b = true;
    const bool c = true;

    PoissonSolver poissonSolver(a,b,c);


    std::array<Real, (NX) * (NY) * (NZ)> F;
    std::array<Real, (NX) * (NY) * (NZ)> sol;

    cout << "Dirichlet test" << endl;
    for (size_t i=0; i < NX; i++){
        for (size_t j=0; j < NY; j++){
            for (size_t k=0; k < NZ; k++){
                F[i * (NY) * (NZ) + j * (NZ) + k] = -3*std::cos(i * DX) *
                                                    std::cos(j * DY) *
                                                    std::sin(k * DZ);

                sol[i * (NY) * (NZ) + j * (NZ) + k] = std::cos(i * DX) *
                                                      std::cos(j * DY) *
                                                      std::sin(k * DZ);
            }   
        }
    }

    fftw_complex* out = fftw_alloc_complex((NX) * (NY) * ((NZ)/2 + 1));

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
/*
    cout << "Neumann test" << endl;
    for (size_t i=0; i < NX+1; i++){
        for (size_t j=0; j < NY+1; j++){
            for (size_t k=0; k < NZ+1; k++){
                F[i * (NY+1) * (NZ+1) + j * (NZ+1) + k] = ((DX/LX) * (DX/LX) + (DY/LY) * (DY/LY) + (DZ/LZ) * (DZ/LZ)) *
                                                          4 * M_PI * M_PI * 
                                                          std::cos(i * DX * 2 * M_PI/LX) *
                                                          std::cos(j * DY * 2 * M_PI/LY) * 
                                                          std::cos(k * DZ * 2 * M_PI/LZ);

                sol[i * (NY+1) * (NZ+1) + j * (NZ+1) + k] = -std::cos(i * DX * 2 * M_PI/LX) *
                                                            std::cos(j * DY * 2 * M_PI/LY) *
                                                            std::cos(k * DZ * 2 * M_PI/LZ);
            }   
        }
    }

    poissonSolver.solveNeumannPoisson(F);

    error = 0;

    for (size_t i=0; i < NX+1; i++){
        for (size_t j=0; j < NY+1; j++){
            for (size_t k=0; k < NZ+1; k++){
                error = error + ((sol[i * (NY+1) * (NZ+1) + j * (NZ+1) + k] - F[i * (NY+1) * (NZ+1) + j * (NZ+1) + k]) *
                                 (sol[i * (NY+1) * (NZ+1) + j * (NZ+1) + k] - F[i * (NY+1) * (NZ+1) + j * (NZ+1) + k]) *
                                 (DX * DY * DZ));
            }   
        }
    }

    std::cout << sqrt(error) << std::endl;
*/
    return 0;
}
