#include "../includes/core.hpp"

int main()
{
    IcoNS icoNS(1.0, 1.0, 1.0, 10, 10, 10, 0.01, 1.0, 100.0, "input.txt", "output.txt");

    icoNS.solve();

    return 0;
}