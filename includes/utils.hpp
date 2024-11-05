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
public:
  double value_x(double x, double y, double z, double t) const
  {
    return std::sin(x * DX) * std::cos(y * DY) * std::sin(z * DZ) * std::sin(t);
  }

  double value_y(double x, double y, double z, double t) const
  {
    return std::cos(x * DX) * std::sin(y * DY) * std::sin(z * DZ) * std::sin(t);
  }

  double value_z(double x, double y, double z, double t) const
  {
    return 2 * std::cos(x * DX) * std::cos(y * DY) * std::cos(z * DZ) * std::sin(t);
  }
};

/**
 * @brief Class containing methods to make controls on our program
 *        Will be unused in a final implementation
 *
 */
class Check
{
private:
  /* data */
public:
  Check(/* args */);

  /**
   * @brief Confronts the state of the grid with the exact solution at a given timestep
   *        How to read: Every cell is rappresented as xyz(GridValue)-(ExactValue)
   *        Printed in 2d squares starting from left face and then moving one step in i direction for each slice of the cube
   *
   * @param grid Whole grid at time t
   * @param exact Exact solution(With DX,DY,DZ correct)
   * @param t Timestep to print at
   * @param d Direction you want to show(U,V,W)
   */
  static void Confront(Grid &grid, ExactSolution &exact, double t, Direction d)
  {
    // U
    if (d == U)
    {
      for (double i = 0; i < NX; i++)
      {
        for (double j = 0; j < NY + 1; j++)
        {
          for (double k = 0; k < NZ + 1; k++)
          {
            std::cout << i << j << k << "(" << grid.u[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] << ")-(" << exact.value_x(i + 0.5, j, k, t) << ") ";
          }
          std::cout << std::endl;
        }
        std::cout << std::endl;
      }
    }
    else if (d == V)
    {
      // V
      for (double i = 0; i < NX + 1; i++)
      {
        for (double j = 0; j < NY; j++)
        {
          for (double k = 0; k < NZ + 1; k++)
          {
            std::cout << i << j << k << "(" << grid.v[i * (NY) * (NZ + 1) + j * (NZ + 1) + k] << ")-(" << exact.value_y(i, j + 0.5, k, t) << ") ";
          }
          std::cout << std::endl;
        }
        std::cout << std::endl;
      }
    }
    else
    {
      // W
      for (double i = 0; i < NX + 1; i++)
      {
        for (double j = 0; j < NY + 1; j++)
        {
          for (double k = 0; k < NZ; k++)
          {
            std::cout << i << j << k << "(" << grid.w[i * (NY + 1) * (NZ) + j * (NZ) + k] << ")-(" << exact.value_z(i, j, k + 0.5, t) << ") ";
          }
          std::cout << std::endl;
        }
        std::cout << std::endl;
      }
    }
  }
};
#endif