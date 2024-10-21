#include <iostream>
#include "core.cpp"

int main()
{
    const double lx = 1000;
    const double ly = 1000;
    const double lz = 1000;
    const unsigned int nx = 100;
    const unsigned int ny = 100;
    const unsigned int nz = 100;
    const double dt = 1;
    const double T = 5.0;
    const double Re = 400.0;
    IcoNS icoNS(lx, ly, lz, nx, ny, nz, dt, T, Re, "input.txt", "output.txt");

    icoNS.solve();

    return 0;
}
