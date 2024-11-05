#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include "utils.hpp"
#include <memory>
#include <vector>

class Boundary
{
    
private:
    //Grid& grid;
    const unsigned int nx;
    const unsigned int ny;
    const unsigned int nz;
    const Real dx,dy,dz;
    //const long Real precision;

    std::vector<std::shared_ptr<BoundaryFunction>> boundary_value_u;
    std::vector<std::shared_ptr<BoundaryFunction>> boundary_value_v;
    std::vector<std::shared_ptr<BoundaryFunction>> boundary_value_w;
    public:
    

    
    Boundary(int nx, int ny, int nz,Real dx_, Real dy_,Real dz_);
    //void update_boundary(long Real t);
    void update_boundary(std::vector<Real>& Yx,std::vector<Real>& Yy,std::vector<Real>& Yz, Real t);

    Real approximate_boundary_u(size_t x, size_t y, size_t z,Real t,size_t face, int side);
    Real approximate_boundary_v(size_t x, size_t y, size_t z,Real t,size_t face, int side);
    Real approximate_boundary_w(size_t x, size_t y, size_t z,Real t,size_t face, int side);

    void addFunction(Direction direction, std::shared_ptr<BoundaryFunction> x);
};
#endif