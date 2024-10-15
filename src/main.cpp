#include <iostream>
#include "core.cpp"

int main()
{
    IcoNS problem(1.0, 1.0, 1.0, 10, 10, 10, 0.01, 1.0, 1000, nullptr, nullptr);
    const std::string input_file = "input.txt";
    problem.preprocessing();
    problem.solve();

    return 0;
}