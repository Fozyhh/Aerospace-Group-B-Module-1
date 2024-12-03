/**
 * @file config.hpp
 * @brief Configuration constants for domain dimensions, discretization, and simulation parameters.
 *
 * This file defines the constants used for setting up the physical domain and the parameters for a
 * numerical simulation, such as the number of cells, space discretization, time step size, and
 * Reynolds number.
 */

#ifndef CONFIG_HPP
#define CONFIG_HPP

#define Real double
#include <cmath>

/// @brief Physical length of the domain in the x-direction.
constexpr Real LX = 2*M_PI;

/// @brief Physical length of the domain in the y-direction.
constexpr Real LY = 2*M_PI;

/// @brief Physical length of the domain in the z-direction.
constexpr Real LZ = 2*M_PI;

/// @brief Number of cells in the x-direction.
constexpr int NX = 5; // 81 cells in each direction gives segmenation fault

/// @brief Number of cells in the y-direction.
constexpr int NY = 5;

/// @brief Number of cells in the z-direction.
constexpr int NZ = 5;

/// @brief Space discretization step size in the x-direction.
constexpr Real DX = LX / NX;

/// @brief Space discretization step size in the y-direction.
constexpr Real DY = LY / NY;

/// @brief Space discretization step size in the z-direction.
constexpr Real DZ = LZ / NZ;

/// @brief Time step size for the simulation.
constexpr Real DT = 0.01;

/// @brief Total time interval for the simulation.
constexpr Real T = 1.0;

/// @brief Reynolds number for the simulation.
constexpr Real RE = 400.0;

#endif // CONFIG_HPP
