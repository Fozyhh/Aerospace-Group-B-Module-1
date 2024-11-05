// class representing the grid of the domain.

#ifndef GRID_HPP
#define GRID_HPP

#include <iostream>
#include <array>
#include "config.hpp"

/**
 * @brief Structure representing the grid of the domain.
 *
 * This class represents the staggered grid used for discretizing the domain. It contains the
 * velocities in the x, y, and z directions, and the pressure field. The grid is staggered in the
 * sense that the velocities are stored at the faces of the cells, while the pressure is stored at
 * grid vertices.
 *
 */
struct Grid
{
    std::array<double, NX *(NY + 1) * (NZ + 1)> u{};
    std::array<double, (NX + 1) * NY *(NZ + 1)> v{};
    std::array<double, (NX + 1) * (NY + 1) * NZ> w{};
    std::array<double, (NX + 1) * (NY + 1) * (NZ + 1)> p{};
};

#endif // GRID_HPP