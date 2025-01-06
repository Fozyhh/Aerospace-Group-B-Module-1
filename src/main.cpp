#include <iostream>
#include "core.hpp"
#include <fstream>
#include <mpi.h>

void writeVTK(const std::string &filename, const double *data, int nx, int ny, int nz, double dx, double dy, double dz)
{
    std::ofstream vtkFile(filename);

    if (!vtkFile.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    // Header
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Navier-Stokes solution\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET STRUCTURED_POINTS\n";
    vtkFile << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n";
    vtkFile << "ORIGIN 0 0 0\n";
    vtkFile << "SPACING " << dx << " " << dy << " " << dz << "\n";

    // Scalar data
    vtkFile << "POINT_DATA " << (nx * ny * nz) << "\n";
    vtkFile << "SCALARS pressure double 1\n";
    vtkFile << "LOOKUP_TABLE default\n";

    for (int k = 0; k < nz; ++k)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                vtkFile << data[i * ny * nz + j * nz + k] << "\n";
            }
        }
    }

    vtkFile.close();
    std::cout << "VTK file written to " << filename << std::endl;
}

void writeVTKVectorField(const std::string &filename, 
                         const double *u, const double *v, const double *w, 
                         int nx, int ny, int nz, double dx, double dy, double dz) 
{
    std::ofstream vtkFile(filename);

    if (!vtkFile.is_open()) 
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    // Header
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Navier-Stokes vector field\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET STRUCTURED_POINTS\n";
    vtkFile << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n";
    vtkFile << "ORIGIN 0 0 0\n";
    vtkFile << "SPACING " << dx << " " << dy << " " << dz << "\n";

    // Vector field data
    vtkFile << "POINT_DATA " << (nx * ny * nz) << "\n";
    vtkFile << "VECTORS velocity double\n";

    for (int k = 0; k < nz; ++k) 
    {
        for (int j = 0; j < ny; ++j) 
        {
            for (int i = 0; i < nx; ++i) 
            {
                int indexx = i * (ny+1) * (nz+1) + j * (nz+1) + k;
                int indexy = i * ny * (nz+1) + j * (nz+1) + k;
                int indexz = i * (ny+1) * nz + j * nz + k;
                vtkFile << u[indexx] << " " << v[indexy] << " " << w[indexz] << "\n";
            }
        }
    }

    vtkFile.close();
    std::cout << "VTK vector field written to " << filename << std::endl;
}


int main(int argc, char* argv[])
{
    
    //Parallel vars
    int rank, size;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double start_time = MPI_Wtime();

    IcoNS icoNS(MPI_COMM_WORLD, argv[1], "output.txt", rank, size);

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

    /* std::array<Real, NX *(NY + 1) * (NZ + 1)> u;
    std::array<Real, (NX + 1) * NY *(NZ + 1)> v;
    std::array<Real, (NX + 1) * (NY + 1) * NZ> w;
    std::array<Real, (NX + 1) * (NY + 1) * (NZ + 1)> p;

    u = icoNS.u_function();
    v = icoNS.v_function();
    w = icoNS.w_function();
    p = icoNS.p_function();

    // Store velocity as a vector field in VTK format
    writeVTKVectorField("velocity.vtk", u.data(), v.data(), w.data(), NX, NY, NZ, DX, DY, DZ);

    // Store pressure as a scalar field in VTK format
    writeVTK("pressure.vtk", p.data(), NX + 1, NY + 1, NZ + 1, DX, DY, DZ);
    writeVTK("u.vtk", u.data(), NX, NY + 1, NZ + 1, DX, DY, DZ);
    writeVTK("v.vtk", v.data(), NX + 1, NY, NZ + 1, DX, DY, DZ);
    writeVTK("w.vtk", w.data(), NX + 1, NY + 1, NZ, DX, DY, DZ);

    fftw_free(icoNS.helper); */
    return 0;
}
