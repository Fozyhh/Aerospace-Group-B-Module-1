#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "core.hpp"
#include "grid.hpp"

class boundary
{
private:
    Grid& grid;
    const unsigned int &nx;
    const unsigned int &ny;
    const unsigned int &nz;
public:
    boundary(Grid& grid_);
    void update_boundary();
    double approximate_boundary_u(size_t x, size_t y, size_t z);
    double boundary_value_u(size_t x, size_t y, size_t z); 
    double boundary_value_v(size_t x, size_t y, size_t z); 
    double boundary_value_w(size_t x, size_t y, size_t z); 
};






#endif