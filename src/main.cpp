#include "../includes/core.hpp"

int main()
{
    IcoNS icoNS(100, 100, 100, 100, 100, 100, 0.01, 10.0, 400.0, "input.txt", "output.txt");

    icoNS.solve();

    return 0;
}