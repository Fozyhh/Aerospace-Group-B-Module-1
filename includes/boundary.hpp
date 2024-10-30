#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include "utils.hpp"
#include "grid.hpp"
#include <memory>
#include <vector>

class Boundary
{
    
private:
    //Grid& grid;
    const unsigned int nx;
    const unsigned int ny;
    const unsigned int nz;
    const double dx,dy,dz;
    const double prec;
    
    std::vector<std::shared_ptr<BoundaryFunction>> boundary_value_u;
    std::vector<std::shared_ptr<BoundaryFunction>> boundary_value_v;
    std::vector<std::shared_ptr<BoundaryFunction>> boundary_value_w;

    
    
    public:

    
    Boundary(int nx, int ny, int nz,double dx_, double dy_,double dz_);
    //void update_boundary(double t);
    void update_boundary(std::vector<double> Yx,std::vector<double> Yy,std::vector<double> Yz, double t);

    double approximate_boundary_u(size_t x, size_t y, size_t z,double t,size_t face);
    double approximate_boundary_v(size_t x, size_t y, size_t z,double t,size_t face);
    double approximate_boundary_w(size_t x, size_t y, size_t z,double t,size_t face);

    void addFunction(size_t direction, std::shared_ptr<BoundaryFunction> x);
};
#endif