/**
 * @file constants.hpp
 * @brief Configuration constants for domain dimensions, discretization, and simulation parameters.
 *
 * This file declares external variables used for setting up the physical domain and the parameters
 * for a numerical simulation. These parameters are defined elsewhere and include domain dimensions,
 * grid resolution, discretization steps, temporal parameters, and flow characteristics.
 */

#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

/**
 * @typedef Real
 * @brief Defines the floating-point precision used throughout the simulation
 *
 * Set to double for double-precision floating-point calculations
 */
#define Real double



extern Real SX,SY,SZ;
/**
 * @brief Physical length of the domain in the x-direction.
 * @details Defines the size of the computational domain along the x-axis in physical units.
 * This value is externally defined and can be modified at runtime.
 */
extern Real LX;

/**
 * @brief Physical length of the domain in the y-direction.
 * @details Defines the size of the computational domain along the y-axis in physical units.
 * This value is externally defined and can be modified at runtime.
 */
extern Real LY;

/**
 * @brief Physical length of the domain in the z-direction.
 * @details Defines the size of the computational domain along the z-axis in physical units.
 * This value is externally defined and can be modified at runtime.
 */
extern Real LZ;

/**
 * @brief Number of cells in the x-direction.
 * @details Specifies the grid resolution along the x-axis. Higher values provide finer resolution
 * but increase computational cost. This value is externally defined and can be modified at runtime.
 */
extern int NX;

/**
 * @brief Number of cells in the y-direction.
 * @details Specifies the grid resolution along the y-axis. Higher values provide finer resolution
 * but increase computational cost. This value is externally defined and can be modified at runtime.
 */
extern int NY;

/**
 * @brief Number of cells in the z-direction.
 * @details Specifies the grid resolution along the z-axis. Higher values provide finer resolution
 * but increase computational cost. This value is externally defined and can be modified at runtime.
 */
extern int NZ;

/**
 * @brief Number of processes in the x-direction for domain decomposition.
 * @details Controls how the domain is split among parallel processes along the x-axis.
 * This value is externally defined and can be modified at runtime.
 */
extern int PX;
/**
 * @brief Number of processes in the y-direction for domain decomposition.
 * @details Controls how the domain is split among parallel processes along the y-axis.
 * This value is externally defined and can be modified at runtime.
 */
extern int PY;
/**
 * @brief Number of processes in the z-direction for domain decomposition.
 * @details Controls how the domain is split among parallel processes along the z-axis.
 * Currently not used as z-direction splitting is not implemented.
 * This value is externally defined and can be modified at runtime.
 */
extern int PZ;

/**
 * @brief Space discretization step size in the x-direction.
 * @details Represents the width of each cell in physical units.
 * Should be calculated as LX/NX. This value is externally defined and can be modified at runtime.
 */
extern Real DX;

/**
 * @brief Space discretization step size in the y-direction.
 * @details Represents the height of each cell in physical units.
 * Should be calculated as LY/NY. This value is externally defined and can be modified at runtime.
 */
extern Real DY;

/**
 * @brief Space discretization step size in the z-direction.
 * @details Represents the depth of each cell in physical units.
 * Should be calculated as LZ/NZ. This value is externally defined and can be modified at runtime.
 */
extern Real DZ;

/**
 * @brief Time step size for the simulation.
 * @details Controls the temporal resolution of the simulation. Smaller values provide more accurate
 * results but increase computational time. Must satisfy CFL condition for stability.
 * This value is externally defined and can be modified at runtime.
 */
extern Real DT;

/**
 * @brief Total number of time steps for the simulation.
 * @details Specifies how long the simulation should run in physical time units.
 * The total simulation time is calculated as Nt * DT.
 * This value is externally defined and can be modified at runtime.
 */
extern Real Nt;

/**
 * @brief Reynolds number for the simulation.
 * @details Dimensionless number that describes the flow regime (laminar vs turbulent).
 * Represents the ratio of inertial forces to viscous forces within the fluid.
 * This value is externally defined and can be modified at runtime.
 */
extern Real RE;

/**
 * @brief Flag for boundary conditions in X direction.
 * @details Specifies the type of boundary conditions to apply at the domain boundaries.
 * false: Periodic boundary conditions
 * true: Dirichlet boundary conditions
 * This value is externally defined and can be modified at runtime.
 */
extern bool BX;

/**
 * @brief Flag for boundary conditions in Z direction.
 * @details Specifies the type of boundary conditions to apply at the domain boundaries.
 * false: Periodic boundary conditions
 * true: Dirichlet boundary conditions
 * This value is externally defined and can be modified at runtime.
 */
extern bool BZ;

/**
 * @brief Flag for boundary conditions in Y direction.
 * @details Specifies the type of boundary conditions to apply at the domain boundaries.
 * false: Periodic boundary conditions
 * true: Dirichlet boundary conditions
 * This value is externally defined and can be modified at runtime.
 */
extern bool BY;

#endif // CONSTANTS_HPP
