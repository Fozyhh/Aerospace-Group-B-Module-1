#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include "utils.hpp"
#include <memory>
#include <vector>

class Boundary
{

private:
    // Grid& grid;

    int lbx,lby,rbx,rby;
    int coords[2];
    std::vector<std::shared_ptr<BoundaryFunction>> boundary_value_u;
    std::vector<std::shared_ptr<BoundaryFunction>> boundary_value_v;
    std::vector<std::shared_ptr<BoundaryFunction>> boundary_value_w;
    int dim_x_x;
    int dim_y_x;
    int dim_x_y;
    int dim_y_y;
    int dim_x_z;
    int dim_y_z;
    int dim_z;
    int dim_z_z;
    int newDimX_x;
    int newDimY_x;
    int newDimX_y;
    int newDimY_y;
    int newDimX_z;
    int newDimY_z;
    int other_dim_x_x;
    int other_dim_y_x;
    int other_dim_x_y;
    int other_dim_y_y;
    int other_dim_x_z;
    int other_dim_y_z;

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

    // void update_boundary(Real t);
    void update_boundary_x(std::vector<Real> &Yx, std::vector<Real> &Yy, std::vector<Real> &Yz, Real t);
    void update_boundary_y(std::vector<Real> &Yx, std::vector<Real> &Yy, std::vector<Real> &Yz, Real t);
    void update_boundary_z(std::vector<Real> &Yx, std::vector<Real> &Yy, std::vector<Real> &Yz, Real t);
    // override
    void update_boundary(std::vector<Real> &Yx, std::vector<Real> &Yy, std::vector<Real> &Yz, Real t);
    Real approximate_boundary_u(size_t x, size_t y, size_t z, Real t, size_t face, int side);
    Real approximate_boundary_v(size_t x, size_t y, size_t z, Real t, size_t face, int side);
    Real approximate_boundary_w(size_t x, size_t y, size_t z, Real t, size_t face, int side);


    void addFunction(Direction direction, std::shared_ptr<BoundaryFunction> x);
    void setBoundaryOffsets(int lbx_, int rbx_, int lby_, int rby_);
    void setCoords(int coords_[2]);
    void setOtherDim(int other_dim_x_x_, int other_dim_y_x_,int other_dim_x_y_, int other_dim_y_y_,int other_dim_x_z_, int other_dim_y_z_);


};
#endif