#include <iostream>
#include "core.hpp"
#include <mpi.h>


int main(int argc, char* argv[])
{
    
    //Parallel vars
    int rank, size;
    MPI_Comm cart_comm;
    int dims[2] = {PX, PY};
    int periods[2] = {1,1};
    
    int coords[2];
    int neighbors[4];

    MPI_Init(&argc, &argv);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Create a Cartesian topology (2D)
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);

    
    MPI_Cart_coords(cart_comm, rank, 2, coords);

    MPI_Cart_shift(cart_comm, 0, 1, &neighbors[0], &neighbors[2]);

    MPI_Cart_shift(cart_comm, 1, 1, &neighbors[1], &neighbors[3]);
    
    IcoNS icoNS("input.txt", "output.txt");

    icoNS.preprocessing();
    icoNS.solve();

    MPI_Finalize();

    return 0;
}
