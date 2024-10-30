#include <iostream>
#include "core.hpp"

int main()
{
    const double lx = 1;
    const double ly = 1;
    const double lz = 1;
    const unsigned int nx = 50;
    const unsigned int ny = 50;
    const unsigned int nz = 50;
    const double dt = 0.001;
    const double T = 1.0;
    const double Re = 400.0;
    IcoNS icoNS(lx, ly, lz, nx, ny, nz, dt, T, Re, "input.txt", "output.txt");
    icoNS.preprocessing();
    icoNS.solve();

    return 0;
}
