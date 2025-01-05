/**
 * @file grid.hpp
 * @brief Defines the staggered grid structure for the fluid simulation domain
 *
 * This file implements a staggered grid arrangement commonly used in computational
 * fluid dynamics. In a staggered grid, velocity components are stored at cell faces
 * while pressure values are stored at cell centers, improving numerical stability
 * and avoiding checker-board pressure patterns.
 */

#ifndef GRID_HPP
#define GRID_HPP

#include <iostream>
#include <vector>


/**
 * @struct Grid
 * @brief Structure representing the staggered grid arrangement for the fluid domain
 *
 * @details
 * The grid uses a staggered arrangement where:
 * - u-velocity components are stored at the x-direction cell faces
 * - v-velocity components are stored at the y-direction cell faces
 * - w-velocity components are stored at the z-direction cell faces
 * - pressure values are stored at cell vertices
 *
 *
 * Array dimensions explanation:
 * - u: NX × (NY+1) × (NZ+1) - stored at x-faces
 * - v: (NX+1) × NY × (NZ+1) - stored at y-faces
 * - w: (NX+1) × (NY+1) × NZ - stored at z-faces
 * - p: (NX+1) × (NY+1) × (NZ+1) - stored at vertices
 *
 */
struct Grid
{
    std::vector<double> u;     ///< x-velocity components
    std::vector<double> v;     ///< y-velocity components
    std::vector<double> w;     ///< z-velocity components
    double *p;                 // z-pencil
    double *px;                // x-pencil
    double *py;                // y-pencil

};

#endif // GRID_HPP
