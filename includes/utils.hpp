#ifndef UTILS_H
#define UTILS_H

#include <functional>
#include <cmath>




class BoundaryFunction{
        public:
            virtual double value(double x, double y, double z,double t) = 0;
};


class FunctionZero : public BoundaryFunction{
    public:
        double value(double x, double y, double z,double t) override {
            return 0;
        }
};

class Dirichlet : public BoundaryFunction{
    public:
    std::function<double(double, double, double ,double)> func;

    Dirichlet(std::function<double(double, double, double ,double)> func_): func(func_){};

    double value(double x, double y, double z,double t) override{
        return func(x,y,z,t);
    }
};

class ExactSolution
  {
    // Place holder implementation for now
  public:
    double value_x(size_t x, size_t y, size_t z, double t) const
    {
      return std::sin(x) * std::cos(y) * std::sin(z) * std::sin(t);
    }

    double value_y(size_t x, size_t y, size_t z, double t) const
    {
      return std::cos(x) * std::sin(y) * std::sin(z) * std::sin(t);
    }

    double value_z(size_t x, size_t y, size_t z, double t) const
    {
      return 2 * std::cos(x) * std::cos(y) * std::cos(z) * std::sin(t);
    }
  };

#endif