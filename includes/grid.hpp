// class representing the grid of the domain.

#ifndef GRID_HPP
#define GRID_HPP

#include <iostream>
#include <array>
#include "config.hpp"

#include <array>

struct Grid
{
    // Define each component as a 3D array with the correct dimensions
    std::array<std::array<std::array<double, NZ + 1>, NY + 1>, NX> u;
    std::array<std::array<std::array<double, NZ + 1>, NY>, NX + 1> v;
    std::array<std::array<std::array<double, NZ>, NY + 1>, NX + 1> w;
    std::array<std::array<std::array<double, NZ + 1>, NY + 1>, NX + 1> p;
};

#endif // GRID_HPP