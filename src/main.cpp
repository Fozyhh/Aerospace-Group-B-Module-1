#include <iostream>
#include "core.hpp"
#include <mpi.h>


int main(int argc, char* argv[])
{

    //Parallel vars
    int rank, size;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double start_time = MPI_Wtime();

    IcoNS icoNS(MPI_COMM_WORLD, "input.txt", "output.txt", rank, size);

    icoNS.preprocessing();
    double solve_start_time = MPI_Wtime();
    icoNS.solve();
    double solve_end_time = MPI_Wtime();
    double end_time = MPI_Wtime();
    if (rank == 0) {
        printf("Total solving time: %f seconds\n", solve_end_time - solve_start_time);
        printf("Total time with init: %f seconds\n", end_time - start_time);
    }

    MPI_Finalize();

    return 0;
}
