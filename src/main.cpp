#include <iostream>
#include "core.hpp"
#include <mpi.h>
#include <iomanip>

void printMemoryRequirements(int rank, int size) {
    if (rank == 0) {
        const std::size_t dim_x_x = NX / PX;
        const std::size_t dim_y_x = (NY + 1) / PY;
        const std::size_t dim_x_y = (NX + 1) / PX;
        const std::size_t dim_y_y = NY / PY;
        const std::size_t dim_x_z = (NX + 1) / PX;
        const std::size_t dim_y_z = (NY + 1) / PY;

        const std::size_t newDimX_x = (dim_x_x + 2);
        const std::size_t newDimY_x = (dim_y_x + 2);
        const std::size_t newDimX_y = (dim_x_y + 2);
        const std::size_t newDimY_y = (dim_y_y + 2);
        const std::size_t newDimX_z = (dim_x_z + 2);
        const std::size_t newDimY_z = (dim_y_z + 2);

        const std::size_t x_size = newDimX_x * newDimY_x * (NZ + 1);
        const std::size_t y_size = newDimX_y * newDimY_y * (NZ + 1);
        const std::size_t z_size = newDimX_z * newDimY_z * NZ;

        const double total_mb = (x_size + y_size + z_size) * sizeof(Real) * 3 / (1024.0 * 1024.0);

        std::cout << std::fixed << std::setprecision(2);
        std::cout << "Memory requirements per process:\n"
                  << "Grid dimensions per process:\n"
                  << "X grid: " << newDimX_x << "x" << newDimY_x << "x" << (NZ+1) << "\n"
                  << "Y grid: " << newDimX_y << "x" << newDimY_y << "x" << (NZ+1) << "\n"
                  << "Z grid: " << newDimX_z << "x" << newDimY_z << "x" << NZ << "\n"
                  << "Total memory per process: " << total_mb << " MB\n"
                  << "Total memory all processes: " << total_mb * size << " MB\n";
    }
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        std::cout << "Starting with " << size << " processes\n";
    }

    // Check process count
    if (size != PX * PY) {
        if (rank == 0) {
            std::cerr << "Error: Number of processes (" << size
                        << ") must equal PX*PY (" << PX << "*" << PY << "="
                        << PX*PY << ")\n";
        }
        MPI_Finalize();
        return 1;
    }

    // Print memory requirements before allocation
    printMemoryRequirements(rank, size);
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "Creating IcoNS instances...\n";
    }

    try {
        IcoNS icoNS(MPI_COMM_WORLD, "input.txt", "output.txt", rank, size);

        if (rank == 0) {
            std::cout << "IcoNS instances created successfully\n";
        }

        std::cout << "Rank " << rank << " completed IcoNS construction\n";
        MPI_Barrier(MPI_COMM_WORLD);

        icoNS.preprocessing();
        std::cout << "Rank " << rank << " completed preprocessing\n";
        MPI_Barrier(MPI_COMM_WORLD);

        icoNS.solve();
        std::cout << "Rank " << rank << " completed solve\n";

    } catch (const std::exception& e) {
        std::cerr << "Rank " << rank << " caught exception: " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }

    MPI_Finalize();
    return 0;
}
