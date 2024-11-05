#include <iostream>
#include "core.hpp"


int main()
{
    const float lx = 1;
    const float ly = 1;
    const float lz = 1;
    const unsigned int nx = 100;
    const unsigned int ny = 100;
    const unsigned int nz = 100;
    const double dt = 0.001;
    const double T = 1.0;
    const float Re = 400.0;
    IcoNS icoNS(lx, ly, lz, nx, ny, nz, dt, T, Re, "input.txt", "output.txt");
    icoNS.preprocessing();
    icoNS.solve();

    return 0;
}
