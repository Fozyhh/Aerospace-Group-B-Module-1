/**
 * @file boundary.hpp
 * @brief Boundary conditions handler
 * @details Manages the application of boundary conditions on all faces of the domain
 */

#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include "utils.hpp"
#include "grid.hpp"
#include <memory>
#include <vector>

/**
 * @class Boundary
 * @brief Handles boundary conditions for the computational domain
 * @details Manages and applies boundary conditions for velocity components (u,v,w)
 *          on all faces of the domain using boundary functions
 */
class Boundary
{
    
private:
    Grid* grid;
    const unsigned int &nx;
    const unsigned int &ny;
    const unsigned int &nz;
    const double dx,dy,dz;
    const double prec;
    
    std::vector<std::shared_ptr<BoundaryFunction>> boundary_value_u;
    std::vector<std::shared_ptr<BoundaryFunction>> boundary_value_v;
    std::vector<std::shared_ptr<BoundaryFunction>> boundary_value_w;

    
    
    public:

    /**
     * @brief Constructs a boundary handler for the given grid
     * @param grid_ Pointer to the grid
     * @param dx_ Grid step in x direction
     * @param dy_ Grid step in y direction
     * @param dz_ Grid step in z direction
     */
    Boundary(Grid* grid_,double dx_, double dy_,double dz_);
    
    /**
     * @brief Updates all boundary values at the given time
     * @param t Current simulation time
     */
    void update_boundary(double t);
    
    /**
     * @brief Approximates boundary value for u-velocity component
     * @param x X-coordinate index
     * @param y Y-coordinate index
     * @param z Z-coordinate index
     * @param t Current time
     * @param face Index of the boundary face
     * @return Approximated boundary value
     */
    double approximate_boundary_u(size_t x, size_t y, size_t z,double t,size_t face);
    
    /**
     * @brief Approximates boundary value for v-velocity component
     * @param x X-coordinate index
     * @param y Y-coordinate index
     * @param z Z-coordinate index
     * @param t Current time
     * @param face Index of the boundary face
     * @return Approximated boundary value
     */
    double approximate_boundary_v(size_t x, size_t y, size_t z,double t,size_t face);
    
    /**
     * @brief Approximates boundary value for w-velocity component
     * @param x X-coordinate index
     * @param y Y-coordinate index
     * @param z Z-coordinate index
     * @param t Current time
     * @param face Index of the boundary face
     * @return Approximated boundary value
     */
    double approximate_boundary_w(size_t x, size_t y, size_t z,double t,size_t face);
    
    /**
     * @brief Adds a boundary function for a specific direction
     * @param direction Direction index for the boundary
     * @param x Shared pointer to the boundary function
     */
    void addFunction(size_t direction, std::shared_ptr<BoundaryFunction> x);
};
#endif
