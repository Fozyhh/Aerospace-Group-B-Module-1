#include <iostream>
#include "core.hpp"

int main()
{
    const double lx = 1000;
    const double ly = 1000;
    const double lz = 1000;
    const unsigned int nx = 64;
    const unsigned int ny = 64;
    const unsigned int nz = 64;
    const double dt = 0.01;
    const double T = 1.0;
    const double Re = 400.0;
    IcoNS icoNS(lx, ly, lz, nx, ny, nz, dt, T, Re, "input.txt", "output.txt");
    icoNS.solve();

    return 0;
}
