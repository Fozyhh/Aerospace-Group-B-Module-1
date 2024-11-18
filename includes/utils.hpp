#ifndef UTILS_H
#define UTILS_H

#include <functional>
#include <cmath>
#include <iostream>
#include <vector>
#include <memory>
#include "grid.hpp"

enum Direction{
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
  // Place holder implementation for now
public:
  ExactSolution(Real dx, Real dy, Real dz) : dx(dx), dy(dy), dz(dz) {}

  Real value_x(Real x, Real y, Real z, Real t) const
  {
    return std::sin(x * dx) * std::cos(y * dy) * std::sin(z * dz) * std::sin(t);
  }

  Real value_y(Real x, Real y, Real z, Real t) const
  {
    return std::cos(x * dx) * std::sin(y * dy) * std::sin(z * dz) * std::sin(t);
  }

  Real value_z(Real x, Real y, Real z, Real t) const
  {
    return 2 * std::cos(x * dx) * std::cos(y * dy) * std::cos(z * dz) * std::sin(t);
  }

private:
  Real dx, dy, dz;
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
   * @param exact Exact solution(With dx,dy,dz correct)
   * @param t Timestep to print at
   * @param d Direction you want to show(U,V,W)
   */
  static void Confront(Grid& grid, ExactSolution& exact, Real t, Direction d){
    //U
    if (d==U){
      for (Real i = 0; i < grid.nx; i++)
      {
          for (Real j = 0; j < grid.ny+1; j++)
          {
              for (Real k = 0; k < grid.nz+1; k++)
              {
                  std::cout << i << j << k << "(" << grid.u[i*(grid.ny+1)*(grid.nz+1) + j*(grid.nz+1) +k] << ")-(" << exact.value_x(i + 0.5,j,k,t) << ") ";
              }
              std::cout << std::endl;
          }
          std::cout << std::endl;
      }
    }else if(d==V){
    //V
      for (Real i = 0; i < grid.nx+1; i++)
      {
          for (Real j = 0; j < grid.ny; j++)
          {
              for (Real k = 0; k < grid.nz+1; k++)
              {
                  std::cout << i << j << k << "(" << grid.v[i*(grid.ny)*(grid.nz+1) + j*(grid.nz+1) +k] << ")-(" << exact.value_y(i,j+0.5,k,t) << ") ";
              }
              std::cout << std::endl;
          }
          std::cout << std::endl;
      }
    }else{
    //W
      for (Real i = 0; i < grid.nx+1; i++)
      {
          for (Real j = 0; j < grid.ny+1; j++)
          {
              for (Real k = 0; k < grid.nz; k++)
              {
                  std::cout << i << j << k << "(" << grid.w[i*(grid.ny+1)*(grid.nz) + j*(grid.nz) +k] << ")-(" << exact.value_z(i,j,k+0.5,t) << ") ";
              }
              std::cout << std::endl;
          }
          std::cout << std::endl;
      }
    }
  }
};
#endif