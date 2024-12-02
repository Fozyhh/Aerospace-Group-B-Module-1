#include <iostream>
#include "core.hpp"
#include <fstream>

void writeVTK(const std::string& filename, const double* data, int nx, int ny, int nz, double dx, double dy, double dz) {
    std::ofstream vtkFile(filename);

    if (!vtkFile.is_open()) {
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
    vtkFile << "SCALARS velocity double 1\n";
    vtkFile << "LOOKUP_TABLE default\n";

    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                vtkFile << data[i * ny * nz + j * nz + k] << "\n";
            }
        }
    }

    vtkFile.close();
    std::cout << "VTK file written to " << filename << std::endl;
}



int main()
{
    IcoNS icoNS("input.txt", "output.txt");

    icoNS.preprocessing();
    icoNS.solve();

    std::array<Real, NX * (NY + 1) * (NZ + 1)> u;
    std::array<Real, (NX + 1) * NY * (NZ + 1)> v;
    std::array<Real, (NX + 1) * (NY + 1) * NZ> w;

    u = icoNS.u_function();
    v = icoNS.v_function();
    w = icoNS.w_function();

    writeVTK("velocity_u.vtk", u.data(), NX, NY + 1, NZ + 1, DX, DY, DZ);
    writeVTK("velocity_v.vtk", v.data(), NX + 1, NY, NZ + 1, DX, DY, DZ);
    writeVTK("velocity_w.vtk", w.data(), NX + 1, NY + 1, NZ, DX, DY, DZ);


    return 0;
}
