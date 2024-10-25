/**
 * @file utils.hpp
 * @brief Utility classes for boundary conditions and exact solutions
 * @details Contains base classes and implementations for boundary conditions
 *          and analytical solutions
 */

#ifndef UTILS_H
#define UTILS_H

#include <functional>
#include <cmath>
#include <iostream>
#include <vector>
#include <memory>

/**
 * @class BoundaryFunction
 * @brief Abstract base class for boundary condition functions
 * @details Defines the interface for all boundary condition implementations
 */
class BoundaryFunction{
    public:
        /**
         * @brief Pure virtual function for computing boundary values
         * @param x X-coordinate
         * @param y Y-coordinate
         * @param z Z-coordinate
         * @param t Time
         * @return Computed boundary value at given position and time
         */
        virtual double value(double x, double y, double z,double t) = 0;
};

/**
 * @class FunctionZero
 * @brief Implements zero boundary condition
 * @details Returns zero for all input coordinates and time
 */
class FunctionZero : public BoundaryFunction{
    public:

        /**
         * @brief Computes zero boundary value
         * @param x X-coordinate (unused)
         * @param y Y-coordinate (unused)
         * @param z Z-coordinate (unused)
         * @param t Time (unused)
         * @return Always returns 0.0
         */
        double value(double x, double y, double z,double t) override {
          return 0;
        }
};

/**
 * @class Dirichlet
 * @brief Implements Dirichlet boundary condition using a function object
 * @details Allows specification of arbitrary function for boundary values
 */
class Dirichlet : public BoundaryFunction{
    public:
        std::function<double(double, double, double ,double)> func;  ///< Function object for boundary values

        /**
         * @brief Constructor taking a function object
         * @param func_ Function object that computes boundary values
         */
        Dirichlet(std::function<double(double, double, double ,double)> func_): func(func_){};

        /**
         * @brief Computes boundary value using stored function
         * @param x X-coordinate
         * @param y Y-coordinate
         * @param z Z-coordinate
         * @param t Time
         * @return Computed boundary value
         */
        double value(double x, double y, double z,double t) override{
          return func(x,y,z,t);
        }
};

/**
 * @class ExactSolution
 * @brief Provides analytical solutions for the problem
 * @details Contains implementations of exact solutions for each velocity component
 */
class ExactSolution{
  public:

      /**
         * @brief Computes exact solution for x-velocity component
         * @param x X-coordinate index
         * @param y Y-coordinate index
         * @param z Z-coordinate index
         * @param t Time
         * @return Exact x-velocity value
         */
      double value_x(size_t x, size_t y, size_t z, double t) const{
        return std::sin(x) * std::cos(y) * std::sin(z) * std::sin(t);
      }

      /**
         * @brief Computes exact solution for y-velocity component
         * @param x X-coordinate index
         * @param y Y-coordinate index
         * @param z Z-coordinate index
         * @param t Time
         * @return Exact y-velocity value
         */
      double value_y(size_t x, size_t y, size_t z, double t) const{
        return std::cos(x) * std::sin(y) * std::sin(z) * std::sin(t);
      }

      /**
         * @brief Computes exact solution for z-velocity component
         * @param x X-coordinate index
         * @param y Y-coordinate index
         * @param z Z-coordinate index
         * @param t Time
         * @return Exact z-velocity value
         */
      double value_z(size_t x, size_t y, size_t z, double t) const{
        return 2 * std::cos(x) * std::cos(y) * std::cos(z) * std::sin(t);
      }
  };


#endif
