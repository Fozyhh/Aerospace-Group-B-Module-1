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
  static void Confront(Grid& grid, ExactSolution& exact, double t, Direction d){
    //U
    if (d==U){
      for (double i = 0; i < grid.nx; i++)
      {
          for (double j = 0; j < grid.ny+1; j++)
          {
              for (double k = 0; k < grid.nz+1; k++)
              {
                  std::cout << i << j << k << "(" << grid.u[i*(grid.ny+1)*(grid.nz+1) + j*(grid.nz+1) +k] << ")-(" << exact.value_x(i + 0.5,j,k,t) << ") ";
              }
              std::cout << std::endl;
          }
          std::cout << std::endl;
      }
    }else if(d==V){
    //V
      for (double i = 0; i < grid.nx+1; i++)
      {
          for (double j = 0; j < grid.ny; j++)
          {
              for (double k = 0; k < grid.nz+1; k++)
              {
                  std::cout << i << j << k << "(" << grid.v[i*(grid.ny)*(grid.nz+1) + j*(grid.nz+1) +k] << ")-(" << exact.value_y(i,j+0.5,k,t) << ") ";
              }
              std::cout << std::endl;
          }
          std::cout << std::endl;
      }
    }else{
    //W
      for (double i = 0; i < grid.nx+1; i++)
      {
          for (double j = 0; j < grid.ny+1; j++)
          {
              for (double k = 0; k < grid.nz; k++)
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