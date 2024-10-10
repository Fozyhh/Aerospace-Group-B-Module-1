// class representing the grid of the domain.

#ifndef GRID_HPP
#define GRID_HPP

#include <iostream>
#include <vector>

struct Grid
{
    int nx, ny, nz;              // number of cells in the x,y,z directions.
    std::vector<double> u, v, w; // velocity values of the grid cells.
    std::vector<double> p;       // pressure values of the grid cells.

    Grid(int nx, int ny, int nz) : nx(nx), ny(ny), nz(nz)
    {
        u.resize(nx * ny * nz);
        v.resize(nx * ny * nz);
        w.resize(nx * ny * nz);
        p.resize((nx + 1) * (ny + 1) * (nz + 1));
    }
};
