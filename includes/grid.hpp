/**
 * @file grid.hpp
 * @brief Grid structure definition
 * @details Defines a 3D grid structure containing velocity and pressure fields
 */

#ifndef GRID_HPP
#define GRID_HPP

#include <iostream>
#include <vector>

/**
 * @class Grid
 * @brief Represents a 3D grid for fluid simulation
 * @details Contains velocity components (u,v,w) and pressure (p) fields stored in 1D vectors
 *          with appropriate indexing for 3D representation
 */
struct Grid
{
    unsigned int nx, ny, nz;     ///< number of cells in the x,y,z directions.
    std::vector<double> u, v, w; ///< velocity values of the grid cells.
    std::vector<double> p;       ///< pressure values of the grid cells.

    /**
     * @brief Constructs a grid with specified dimensions
     * @param nx Number of cells in x direction
     * @param ny Number of cells in y direction
     * @param nz Number of cells in z direction
     */
    Grid(int nx, int ny, int nz) : nx(nx), ny(ny), nz(nz)
    {
        u.resize(nx * (ny+1) * (nz+1));
        v.resize((nx+1) * ny * (nz+1));
        w.resize((nx+1) * (ny+1) * nz);
        p.resize((nx + 1) * (ny + 1) * (nz + 1));
    }
};

#endif // GRID_HPP
