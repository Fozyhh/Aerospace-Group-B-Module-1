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

    IcoNS icoNS(MPI_COMM_WORLD, argv[1], rank, size);

    icoNS.solve();

    double end_time = MPI_Wtime();

    if (rank == 0) {
        std::cout << "Total time with init and output: " << end_time - start_time << " seconds" << std::endl;
    }


    MPI_Finalize();

    return 0;
}
