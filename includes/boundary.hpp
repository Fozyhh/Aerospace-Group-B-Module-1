#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include "utils.hpp"
#include <memory>
#include <vector>

/**
 * @class Boundary
 * @brief Class handling boundary conditions for a 3D fluid simulation
 *
 * This class manages boundary conditions for velocity components (u, v, w) in a
 * three-dimensional fluid simulation. It provides functionality for updating boundary
 * values, calculating approximate boundary values, and managing boundary functions.
 */
class Boundary
{
private:
    /// @brief Left and right boundary offsets in x and y directions
    int lbx,lby,rbx,rby;

    /// @brief Process coordinates in the 2D grid
    int coords[2];

    /// @brief Vector of boundary functions for u-velocity component
    std::vector<std::shared_ptr<BoundaryFunction>> boundary_value_u;
    /// @brief Vector of boundary functions for v-velocity component
    std::vector<std::shared_ptr<BoundaryFunction>> boundary_value_v;
    /// @brief Vector of boundary functions for w-velocity component
    std::vector<std::shared_ptr<BoundaryFunction>> boundary_value_w;

    /// @brief Dimension parameters for x-direction mesh
    int dim_x_x, dim_y_x;
    /// @brief Dimension parameters for y-direction mesh
    int dim_x_y, dim_y_y, dim_y_z;
    /// @brief Dimension parameters for z-direction mesh
    int dim_z, dim_x_z, dim_z_z;

    /// @brief New dimension parameters for x-direction mesh
    int newDimX_x, newDimY_x;
    /// @brief New dimension parameters for y-direction mesh
    int newDimX_y, newDimY_y;
    /// @brief New dimension parameters for z-direction mesh
    int newDimX_z, newDimY_z;

    int zSize[3];

    /// @brief Additional dimension parameters for x-direction mesh
    int other_dim_x_x, other_dim_y_x;
    /// @brief Additional dimension parameters for y-direction mesh
    int other_dim_x_y, other_dim_y_y;
    /// @brief Additional dimension parameters for z-direction mesh
    int other_dim_x_z, other_dim_y_z;

public:
    /**
    * @brief Initialize boundary dimensions and parameters
    * @param dim_x_x_ X dimension for x-direction mesh
    * @param dim_y_x_ Y dimension for x-direction mesh
    * @param dim_x_y_ X dimension for y-direction mesh
    * @param dim_y_y_ Y dimension for y-direction mesh
    * @param dim_x_z_ X dimension for z-direction mesh
    * @param dim_y_z_ Y dimension for z-direction mesh
    * @param dim_z_ Z dimension
    * @param dim_z_z_ Additional Z dimension parameter
    * @param newDimX_x_ New X dimension for x-direction mesh
    * @param newDimY_x_ New Y dimension for x-direction mesh
    * @param newDimX_y_ New X dimension for y-direction mesh
    * @param newDimY_y_ New Y dimension for y-direction mesh
    * @param newDimX_z_ New X dimension for z-direction mesh
    * @param newDimY_z_ New Y dimension for z-direction mesh
    */
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

    //TODO: check, Y2 to change in vector most likely
    void divergence(std::vector<Real> &Yx, std::vector<Real> &Yy, std::vector<Real> &Yz, std::vector<Real> &Y2_p/*std::array<Real, (NX + 1) * (NY + 1) * (NZ + 1)> &Y2_p*/, Real t, Real c);

    /**
     * @brief Adds a boundary function for a specified velocity direction
     * @param direction Velocity direction (U=length, V=width, W=height)
     * @param x Shared pointer to the boundary function
     */
    void addFunction(Direction direction, std::shared_ptr<BoundaryFunction> x);

    /**
     * @brief Sets boundary offsets
     * @param lbx_ Left boundary offset in x-direction
     * @param rbx_ Right boundary offset in x-direction
     * @param lby_ Left boundary offset in y-direction
     * @param rby_ Right boundary offset in y-direction
    */
    void setBoundaryOffsets(int lbx_, int rbx_, int lby_, int rby_);

    /**
     * @brief Sets process coordinates in the 2D grid
     * @param coords_ Array containing process coordinates
    */
    void setCoords(int coords_[2]);

    /**
     * @brief Sets additional dimension parameters
     * @param other_dim_x_x_ Additional X dimension for x-direction mesh
     * @param other_dim_y_x_ Additional Y dimension for x-direction mesh
     * @param other_dim_x_y_ Additional X dimension for y-direction mesh
     * @param other_dim_y_y_ Additional Y dimension for y-direction mesh
     * @param other_dim_x_z_ Additional X dimension for z-direction mesh
     * @param other_dim_y_z_ Additional Y dimension for z-direction mesh
    */
    void setOtherDim(int other_dim_x_x_, int other_dim_y_x_,int other_dim_x_y_, int other_dim_y_y_,int other_dim_x_z_, int other_dim_y_z_);


};
#endif
