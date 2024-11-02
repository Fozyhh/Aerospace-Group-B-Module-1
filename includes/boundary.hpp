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
    // void update_boundary(double t);
    void update_boundary(std::vector<double> &Yx, std::vector<double> &Yy, std::vector<double> &Yz, double t);
    // override
    void update_boundary(std::array<double, NX *(NY + 1) * (NZ + 1)> &Yx, std::array<double, (NX + 1) * NY *(NZ + 1)> &Yy, std::array<double, (NX + 1) * (NY + 1) * NZ> &Yz, double t);

    double approximate_boundary_u(size_t x, size_t y, size_t z, double t, size_t face, int side);
    double approximate_boundary_v(size_t x, size_t y, size_t z, double t, size_t face, int side);
    double approximate_boundary_w(size_t x, size_t y, size_t z, double t, size_t face, int side);

    void addFunction(Direction direction, std::shared_ptr<BoundaryFunction> x);
};
#endif