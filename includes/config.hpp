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

/// @brief Physical length of the domain in the x-direction.
constexpr Real LX = 1.0;

/// @brief Physical length of the domain in the y-direction.
constexpr Real LY = 1.0;

/// @brief Physical length of the domain in the z-direction.
constexpr Real LZ = 1.0;

/// @brief Number of cells in the x-direction.
constexpr int NX = 50;

/// @brief Number of cells in the y-direction.
constexpr int NY = 50;

/// @brief Number of cells in the z-direction.
constexpr int NZ = 50;

constexpr int PX=2;
constexpr int PY=1;
constexpr int PZ=1;//Probably useless as we dont want to split z

constexpr int dim_x_x = NX / PX;
constexpr int dim_y_x = (NY + 1) / PY;
constexpr int dim_x_y = (NX + 1) / PX;
constexpr int dim_y_y = NY / PY;
constexpr int dim_x_z = (NX + 1) / PX;
constexpr int dim_y_z = (NY + 1) / PY;
constexpr int dim_z = NZ + 1;
constexpr int dim_z_z = NZ;
constexpr int newDimX_x = (dim_x_x+2);
constexpr int newDimY_x = (dim_y_x+2);
constexpr int newDimX_y = (dim_x_y+2);
constexpr int newDimY_y = (dim_y_y+2);
constexpr int newDimX_z = (dim_x_z+2);
constexpr int newDimY_z = (dim_y_z+2);



/// @brief Space discretization step size in the x-direction.
constexpr Real DX = LX / NX;

/// @brief Space discretization step size in the y-direction.
constexpr Real DY = LY / NY;

/// @brief Space discretization step size in the z-direction.
constexpr Real DZ = LZ / NZ;

/// @brief Time step size for the simulation.
constexpr Real DT = 0.001;

/// @brief Total time interval for the simulation.
constexpr Real T = 0.5;

/// @brief Reynolds number for the simulation.
constexpr Real RE = 1000.0;

#endif // CONFIG_HPP
