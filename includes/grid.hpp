// class representing the grid of the domain.

#ifndef GRID_HPP
#define GRID_HPP

#include <iostream>
#include <array>
#include "config.hpp"

struct Grid
{
    std::array<double, NX *(NY + 1) * (NZ + 1)> u{};
    std::array<double, (NX + 1) * NY *(NZ + 1)> v{};
    std::array<double, (NX + 1) * (NY + 1) * NZ> w{};
    std::array<double, (NX + 1) * (NY + 1) * (NZ + 1)> p{};
};

#endif // GRID_HPP