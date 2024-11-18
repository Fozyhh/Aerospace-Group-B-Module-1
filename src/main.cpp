#include <iostream>
#include "core.hpp"


int main()
{
    const float lx = 1;
    const float ly = 1;
    const float lz = 1;
    const unsigned int nx = 4;
    const unsigned int ny = 4;
    const unsigned int nz = 4;
    const double dt = 0.01;
    const double T = 1.0;
    const float Re = 400.0;
    IcoNS icoNS(lx, ly, lz, nx, ny, nz, dt, T, Re, "input.txt", "output.txt","error_log.txt");
    icoNS.preprocessing();
    icoNS.solve();

    return 0;
}
