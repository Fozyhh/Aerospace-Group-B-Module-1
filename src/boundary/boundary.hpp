#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "core.hpp"
#include "grid.hpp"
#include <vector>
#include <functional>

class boundary
{
private:
    Grid& grid;
    const unsigned int &nx;
    const unsigned int &ny;
    const unsigned int &nz;
    const double dx,dy,dz;
    const double prec;
    std::vector<Function> boundary_value_u;
    std::vector<Function> boundary_value_v;
    std::vector<Function> boundary_value_w;
    
public:
    boundary(Grid& grid_,double dx_, double dy_,double dz_);
    void update_boundary(double t);

    double boundary::approximate_boundary_u(size_t x, size_t y, size_t z,double t,size_t face);
    double boundary::approximate_boundary_v(size_t x, size_t y, size_t z,double t,size_t face);
    double boundary::approximate_boundary_w(size_t x, size_t y, size_t z,double t,size_t face);
    //double boundary_value_u(size_t x, size_t y, size_t z); 
    //double boundary_value_v(size_t x, size_t y, size_t z); 
    //double boundary_value_w(size_t x, size_t y, size_t z);

    void addFunction(size_t direction, Function& x);

};

class Function
{
public:
    virtual double value(double x, double y, double z,double t) = 0;
};


class FunctionZero : public Function{
public:
    double value(double x, double y, double z,double t) override {
        return 0;
    }
};

class FunctionX : public Function{
    public:
    std::function<double(double, double, double ,double)> func;

    FunctionX(std::function<double(double, double, double ,double)> func_): func(func_){};

    double value(double x, double y, double z,double t) override{
        return func(x,y,z,t);
    }
}







#endif