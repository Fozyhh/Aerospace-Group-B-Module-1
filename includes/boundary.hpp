#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include "utils.hpp"
#include <memory>
#include <vector>

/**
 * @brief //TODO
 */
class Boundary
{
private:

    int lbx,lby,rbx,rby;
    int coords[2];

    std::vector<std::shared_ptr<BoundaryFunction>> boundary_value_u;
    std::vector<std::shared_ptr<BoundaryFunction>> boundary_value_v;
    std::vector<std::shared_ptr<BoundaryFunction>> boundary_value_w;

    int dim_x_x, dim_y_x;
    int dim_x_y, dim_y_y, dim_y_z;
    int dim_z, dim_x_z, dim_z_z;

    int newDimX_x, newDimY_x; 
    int newDimX_y, newDimY_y;
    int newDimX_z, newDimY_z;
    
    int other_dim_x_x, other_dim_y_x;
    int other_dim_x_y, other_dim_y_y;
    int other_dim_x_z, other_dim_y_z;

public:
    void initializeBoundary(
        int dim_x_x_, int dim_y_x_, int dim_x_y_, int dim_y_y_,
        int dim_x_z_, int dim_y_z_, int dim_z_, int dim_z_z_,
        int newDimX_x_, int newDimY_x_, int newDimX_y_, int newDimY_y_,
        int newDimX_z_, int newDimY_z_)
    {
        dim_x_x = dim_x_x_;
        dim_y_x = dim_y_x_;
        dim_x_y = dim_x_y_;
        dim_y_y = dim_y_y_;
        dim_x_z = dim_x_z_;
        dim_y_z = dim_y_z_;
        dim_z = dim_z_;
        dim_z_z = dim_z_z_;
        newDimX_x = newDimX_x_;
        newDimY_x = newDimY_x_;
        newDimX_y = newDimX_y_;
        newDimY_y = newDimY_y_;
        newDimX_z = newDimX_z_;
        newDimY_z = newDimY_z_; 
    }
 
    /**
     * @brief The method is called by the program multiple during the time step, 
     *        in order to update the values of the boundaries at each
     *        requested time t, calculating the approximated ones too.
     *
     * @param Yx Boundary x velocities or the Y intermediate function related to the x direction.
     * @param Yy Boundary y velocities or the Y intermediate function related to the y direction.
     * @param Yz Boundary z velocities or the Y intermediate function related to the z direction.
     * @param t Time of the time discretization we are considering.
     */
    void update_boundary(std::vector<Real> &Yx, std::vector<Real> &Yy, std::vector<Real> &Yz, Real t);

    /**
     * @brief Calculate the approximate value of the x velocity in a given point.
     *
     * @param x,y,z Coordinates of the position in the 3D mesh.
     * @param t Time of the time discretization we are considering.
     * @param face Face of the mesh we are considering. In order: 0 (LEFT), 1 (RIGHT), 2 (FRONT), 3 (BACK), 4 (LOWER), 5 (UPPER)
     * @param side Int used to calculate the correct approximate value. It's 1 for faces 0, 2 and 4, -1 for the others.
     *
     * @return the approximate value.
     */
    Real approximate_boundary_u(int x, int y, int z, Real t, int face, int side);

    
    /**
     * @brief Calculate the approximate value of the y velocity in a given point.
     *
     * @param x,y,z Coordinates of the position in the 3D mesh.
     * @param t Time of the time discretization we are considering.
     * @param face Face of the mesh we are considering. In order: 0 (LEFT), 1 (RIGHT), 2 (FRONT), 3 (BACK), 4 (LOWER), 5 (UPPER)
     * @param side Int used to calculate the correct approximate value. It's 1 for faces 0, 2 and 4, -1 for the others.
     *
     * @return the approximate value.
     */
    Real approximate_boundary_v(int x, int y, int z, Real t, int face, int side);

    /**
     * @brief Calculate the approximate value of the z velocity in a given point.
     *
     * @param x,y,z Coordinates of the position in the 3D mesh.
     * @param t Time of the time discretization we are considering.
     * @param face Face of the mesh we are considering. In order: 0 (LEFT), 1 (RIGHT), 2 (FRONT), 3 (BACK), 4 (LOWER), 5 (UPPER)
     * @param side Int used to calculate the correct approximate value. It's 1 for faces 0, 2 and 4, -1 for the others.
     *
     * @return the approximate value.
     */
    Real approximate_boundary_w(int x, int y, int z, Real t, int face, int side);

    /**
     * @brief Add the given function to the selected direction.
     *
     * @param direction Direction U (length), V (width) or W (height) of the boundary.
     * @param x Function to assign to the boundary
     */
    void addFunction(Direction direction, std::shared_ptr<BoundaryFunction> x);

    
    void setBoundaryOffsets(int lbx_, int rbx_, int lby_, int rby_);
    void setCoords(int coords_[2]);
    void setOtherDim(int other_dim_x_x_, int other_dim_y_x_,int other_dim_x_y_, int other_dim_y_y_,int other_dim_x_z_, int other_dim_y_z_);


};
#endif