#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include "utils.hpp"
#include <memory>
#include <vector>

class Boundary
{

private:
    // Grid& grid;


    std::vector<std::shared_ptr<BoundaryFunction>> boundary_value_u;
    std::vector<std::shared_ptr<BoundaryFunction>> boundary_value_v;
    std::vector<std::shared_ptr<BoundaryFunction>> boundary_value_w;


public:
    // void update_boundary(Real t);
    void update_boundary(std::vector<Real> &Yx, std::vector<Real> &Yy, std::vector<Real> &Yz, std::vector<Real> &P, Real t);
    // override
    void update_boundary(std::array<Real, NX *(NY + 1) * (NZ + 1)> &Yx, std::array<Real, (NX + 1) * NY *(NZ + 1)> &Yy, std::array<Real, (NX + 1) * (NY + 1) * NZ> &Yz, std::array<Real, (NX + 1) * (NY + 1) * (NZ + 1)> &P, Real t);

    Real approximate_boundary_u(int x, int y, int z, Real t, int face, int side);
    Real approximate_boundary_v(int x, int y, int z, Real t, int face, int side);
    Real approximate_boundary_w(int x, int y, int z, Real t, int face, int side);


    void addFunction(Direction direction, std::shared_ptr<BoundaryFunction> x);
};
#endif