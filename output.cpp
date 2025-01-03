#include <iostream>
#include <array>
#include <mpi.h>
#include <vector>

#define NX 4
#define NY 4
#define NZ 4

#define PX 2
#define PY 1
#define PZ 1

int main(int argc, char **argv)
{

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Create a Cartesian topology (2D)
    MPI_Comm cart_comm;
    int dims[2] = {PX, PY};
    int periods[2] = {0,0};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);

    int dim_x = NX / PX;
    int dim_y = (NY + 1) / PY;
    int dim_z = NZ + 1;
    int coords[2];
    MPI_Cart_coords(cart_comm, rank, 2, coords);
    int neighbors[4];

    MPI_Cart_shift(cart_comm, 0, 1, &neighbors[0], &neighbors[2]);

    MPI_Cart_shift(cart_comm, 1, 1, &neighbors[1], &neighbors[3]);

    if (NX % PX != 0 && coords[0] == dim_x - 1)
        dim_x++;
    if ((NY + 1) % PY != 0 && coords[1] == 1/*dim_y - 1*/){
        dim_y++;
    }
        
    int res_x = 0;
    int res_y = 0;

    for (int i = 0; i < 4; i++)
    {
        if (neighbors[i] != -2)
        {
            if (i == 0 || i == 2)
            {
                res_x++;
            }
            else
            {
                res_y++;
            }
        }
    }
    int newDimX = dim_x + res_x;
    int newDimY = dim_y + res_y;
    std::vector<int> grid_loc(newDimX * newDimY * dim_z);


    // Write output
    const char *filename = "solution.vtk";

    // Open file
    MPI_File fh;
    MPI_File_open(cart_comm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    MPI_Offset headerOffset = 0;
    if (rank == 0) {
        // Write header
        std::ostringstream header;
        header << "# vtk DataFile Version 3.0\n"
               << "Aerospace Module 1 Output\n"
               << "BINARY\n"
               << "DATASET STRUCTURED_POINTS\n"
               << "DIMENSIONS " << Nx << " " << Ny << " 1\n"
               << "ORIGIN " << 0 << " " << 0 << " " << 0 << "\n"
               << "SPACING " << dx << " " << dy << " " << dz << "\n"
               << "POINT_DATA " << Nx * Ny << "\n";
        std::string header_str = header.str();
        MPI_File_write_at(fh, 0, header_str.c_str(), header_str.size(), MPI_CHAR, MPI_STATUS_IGNORE);
        headerOffset = header_str.size();
    }

    // Write data
    MPI_Offset dataOffset = headerOffset + rank * newDimX * newDimY * dim_z * sizeof(int);






    MPI_File_close(&fh);

    MPI_Finalize();
    return 0;
}
