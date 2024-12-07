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


public:
    // void update_boundary(Real t);
    void update_boundary(std::vector<Real> &Yx, std::vector<Real> &Yy, std::vector<Real> &Yz, Real t);
    // override
    void update_boundary(std::array<Real,newDimX_x*newDimY_x*dim_z> &Yx, std::array<Real, newDimX_y*newDimY_y*dim_z> &Yy, std::array<Real, newDimX_z*newDimY_z*dim_z_z> &Yz, Real t);
    Real approximate_boundary_u(size_t x, size_t y, size_t z, Real t, size_t face, int side);
    Real approximate_boundary_v(size_t x, size_t y, size_t z, Real t, size_t face, int side);
    Real approximate_boundary_w(size_t x, size_t y, size_t z, Real t, size_t face, int side);


    void addFunction(Direction direction, std::shared_ptr<BoundaryFunction> x);
    void setBoundaryOffsets(int lbx_, int rbx_, int lby_, int rby_);
    void setCoords(int coords[2]);
};
#endif