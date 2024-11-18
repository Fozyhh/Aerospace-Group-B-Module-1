#include <iostream>
#include "core.hpp"


int main()
{
    const float lx = 1;
    const float ly = 1;
    const float lz = 1;
    const double T = 0.5;
    const float Re = 400.0;

    IcoNS icoNS10(lx, ly, lz, 10, 10, 10, 1/100.0, T, Re, "input.txt", "output.txt","conv/10.txt");
    icoNS10.preprocessing();
    icoNS10.solve();

    IcoNS icoNS20(lx, ly, lz, 20, 20, 20, 1/200.0, T, Re, "input.txt", "output.txt","conv/20.txt");
    icoNS20.preprocessing();
    icoNS20.solve();

    IcoNS icoNS50(lx, ly, lz, 50, 50, 50, 1/500.0, T, Re, "input.txt", "output.txt","conv/50.txt");
    icoNS50.preprocessing();
    icoNS50.solve();

    IcoNS icoNS100(lx, ly, lz, 100, 100, 100, 1/1000.0, T, Re, "input.txt", "output.txt","conv/z100.txt");
    icoNS100.preprocessing();
    icoNS100.solve();
    return 0;
}
