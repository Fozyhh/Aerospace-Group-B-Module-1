#ifndef UTILS_H
#define UTILS_H

#include <functional>
#include <cmath>
#include <iostream>
#include <vector>
#include <memory>

class BoundaryFunction
{
public:
  virtual double value(double x, double y, double z, double t) = 0;
};

class FunctionZero : public BoundaryFunction
{
public:
  double value(double /*x*/, double /*y*/, double /*z*/, double /*t*/) override
  {
    return 0;
  }
};

class Dirichlet : public BoundaryFunction
{
public:
  std::function<double(double, double, double, double)> func;

  Dirichlet(std::function<double(double, double, double, double)> func_) : func(func_) {};

  double value(double x, double y, double z, double t) override
  {
    return func(x, y, z, t);
  }
};

class ExactSolution
{
  // Place holder implementation for now
public:
  ExactSolution(double dx, double dy, double dz) : dx(dx), dy(dy), dz(dz) {}

  double value_x(double x, double y, double z, double t) const
  {
    return std::sin(x * dx) * std::cos(y * dy) * std::sin(z * dz) * std::sin(t);
  }

  double value_y(double x, double y, double z, double t) const
  {
    return std::cos(x * dx) * std::sin(y * dy) * std::sin(z * dz) * std::sin(t);
  }

  double value_z(double x, double y, double z, double t) const
  {
    return 2 * std::cos(x * dx) * std::cos(y * dy) * std::cos(z * dz) * std::sin(t);
  }

private:
  double dx, dy, dz;
};

#endif