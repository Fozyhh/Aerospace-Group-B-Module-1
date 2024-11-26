#include <iostream>
#include "core.hpp"


int main()
{
    const float lx = 1;
    const float ly = 1;
    const float lz = 1;
    const double T = 0.5;
    const float Re = 400.0;

    IcoNS icoNS10("input.txt", "output.txt");
    icoNS10.preprocessing();
    icoNS10.solve();

    IcoNS icoNS20("input.txt", "output.txt");
    icoNS20.preprocessing();
    icoNS20.solve();

    IcoNS icoNS50("input.txt", "output.txt");
    icoNS50.preprocessing();
    icoNS50.solve();

    IcoNS icoNS100("input.txt", "output.txt");
    icoNS100.preprocessing();
    icoNS100.solve();
    return 0;
}
