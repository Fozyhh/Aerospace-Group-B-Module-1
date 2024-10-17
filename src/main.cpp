#include <iostream>
#include "core.cpp"

int main()
{
    IcoNS problem(1.0, 1.0, 1.0, 50, 50, 50, 0.001, 1.0, 1000, "", "");
    const std::string input_file = "input.txt";
    problem.preprocessing();
    problem.solve();

    return 0;
}