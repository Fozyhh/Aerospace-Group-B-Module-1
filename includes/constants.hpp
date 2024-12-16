/**
 * @file constants.hpp
 * @brief Configuration constants for domain dimensions, discretization, and simulation parameters.
 *
 * This file defines the constants used for setting up the physical domain and the parameters for a
 * numerical simulation, such as the number of cells, space discretization, time step size, and
 * Reynolds number.
 */

#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#define Real double

/// @brief Physical length of the domain in the x-direction.
extern Real LX;

/// @brief Physical length of the domain in the y-direction.
extern Real LY;

/// @brief Physical length of the domain in the z-direction.
extern Real LZ;

/// @brief Number of cells in the x-direction.
extern int NX;

/// @brief Number of cells in the y-direction.
extern int NY;

/// @brief Number of cells in the z-direction.
extern int NZ;

extern int PX;
extern int PY;
extern int PZ; // Probably useless as we dont want to split z

/// @brief Space discretization step size in the x-direction.
extern Real DX;

/// @brief Space discretization step size in the y-direction.
extern Real DY;

/// @brief Space discretization step size in the z-direction.
extern Real DZ;

/// @brief Time step size for the simulation.
extern Real DT;

/// @brief Total time interval for the simulation.
extern Real T;

/// @brief Reynolds number for the simulation.
extern Real RE;

#endif // CONSTANTS_HPP
