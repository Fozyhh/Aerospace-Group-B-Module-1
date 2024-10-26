#include <iostream>
#include "core.hpp"

int main()
{
    const double lx = 1000;
    const double ly = 1000;
    const double lz = 1000;
    const double dt = 0.1;
    const double T = 5.0;
    const double Re = 400.0;
    IcoNS icoNS(lx, ly, lz, dt, T, Re, "input.txt", "output.txt");

    icoNS.solve();

    return 0;
}
