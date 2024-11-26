#include <iostream>
#include "core.hpp"


int main()
{
    IcoNS icoNS("input.txt", "output.txt");

    icoNS.preprocessing();
    icoNS.solve();

    return 0;
}
