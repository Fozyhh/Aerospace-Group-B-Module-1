#ifndef UTILS_H
#define UTILS_H

#include <functional>
#include <cmath>
#include <iostream>
#include <vector>
#include <memory>
#include "grid.hpp"

enum Direction
{
  U,
  V,
  W
};

class BoundaryFunction
{
public:
  virtual Real value(Real x, Real y, Real z, Real t) = 0;
};

class FunctionZero : public BoundaryFunction
{
public:
  Real value(Real /*x*/, Real /*y*/, Real /*z*/, Real /*t*/) override
  {
    return 0;
  }
};

class Dirichlet : public BoundaryFunction
{
public:
  std::function<Real(Real, Real, Real, Real)> func;

  Dirichlet(std::function<Real(Real, Real, Real, Real)> func_) : func(func_) {};

  Real value(Real x, Real y, Real z, Real t) override
  {
    return func(x, y, z, t);
  }
};

class ExactSolution
{
public:

  Real value_x(Real x, Real y, Real z, Real t) const
  {
    return std::sin(x * DX) * std::cos(y * DY) * std::sin(z * DZ) * std::sin(t);
  }

  Real value_y(Real x, Real y, Real z, Real t) const
  {
    return std::cos(x * DX) * std::sin(y * DY) * std::sin(z * DZ) * std::sin(t);
  }

  Real value_z(Real x, Real y, Real z, Real t) const
  {
    return 2 * std::cos(x * DX) * std::cos(y * DY) * std::cos(z * DZ) * std::sin(t);
  }
};
#endif