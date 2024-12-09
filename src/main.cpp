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
    
    IcoNS icoNS(MPI_COMM_WORLD, "input.txt", "output.txt", rank, size);

    icoNS.preprocessing();
    icoNS.solve();

    MPI_Finalize();

    return 0;
}
