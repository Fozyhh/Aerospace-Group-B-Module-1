/**
 * @file utils.hpp
 * @brief Utility classes and functions for boundary conditions and exact solutions
 *
 * This file provides utility classes for handling boundary conditions and exact solutions
 * in the fluid simulation. It includes enumerations for direction, abstract and concrete
 * boundary condition classes, and exact solution calculations.
 */

#ifndef UTILS_H
#define UTILS_H

#include <functional>
#include <cmath>
#include <iostream>
#include <vector>
#include <memory>
#include <string>
#include <sstream>
#include "grid.hpp"
#include "constants.hpp"
#include <iomanip>

/**
 * @enum Direction
 * @brief Enumeration for velocity components direction
 *
 * @var Direction::U X-direction velocity component
 * @var Direction::V Y-direction velocity component
 * @var Direction::W Z-direction velocity component
 */
enum Direction
{
    U,  ///< X-direction
    V,  ///< Y-direction
    W   ///< Z-direction
};

enum Faces
{
  LEFT,
  RIGHT,
  FRONT,
  BACK,
  LOWER,
  UPPER
};




/**
 * @class BoundaryFunction
 * @brief Abstract base class for boundary condition functions
 *
 * Defines the interface for boundary condition functions. All specific boundary
 * condition types must inherit from this class and implement the value method.
 */
class BoundaryFunction
{
public:
    /**
    * @brief Pure virtual function to calculate boundary value
    * @param x X-coordinate
    * @param y Y-coordinate
    * @param z Z-coordinate
    * @param t Time
    * @return Boundary value at the specified position and time
    */
    virtual Real value(Real x, Real y, Real z, Real t) = 0;
};

/**
 * @class FunctionZero
 * @brief Implements zero boundary condition
 *
 * Derived class from BoundaryFunction that returns zero regardless of position or time.
 * Useful for no-slip or zero-flux boundary conditions.
 */
class FunctionZero : public BoundaryFunction
{
public:
  /**
   * @brief Implements zero boundary condition
   * @param x X-coordinate (unused)
   * @param y Y-coordinate (unused)
   * @param z Z-coordinate (unused)
   * @param t Time (unused)
   * @return Always returns 0
   */
  Real value(Real /*x*/, Real /*y*/, Real /*z*/, Real /*t*/) override
  {
    return 0;
  }
};

/**
 * @class Dirichlet
 * @brief Implements Dirichlet boundary condition
 *
 * Allows specification of arbitrary function for Dirichlet boundary conditions.
 * The boundary value is determined by the provided function.
 */
class Dirichlet : public BoundaryFunction
{
public:
  /// Function object that defines the boundary value
  std::function<Real(Real, Real, Real, Real)> func;

  /**
    * @brief Constructor taking a function object
    * @param func_ Function object defining boundary values
    */
  Dirichlet(std::function<Real(Real, Real, Real, Real)> func_) : func(func_) {};

  /**
    * @brief Calculates boundary value using the stored function
    * @param x X-coordinate
    * @param y Y-coordinate
    * @param z Z-coordinate
    * @param t Time
    * @return Boundary value calculated by the stored function
    */
  Real value(Real x, Real y, Real z, Real t) override
  {
    return func(x, y, z, t);
  }
};

/**
 * @class ExactSolution
 * @brief Provides exact solutions for velocity components
 *
 * Contains analytical solutions for velocity components in x, y, and z directions.
 * Used for validation and error analysis of numerical solutions.
 */
class ExactSolution
{
public:

  /**
    * @brief Calculates exact solution for x-velocity component
    * @param x X-coordinate
    * @param y Y-coordinate
    * @param z Z-coordinate
    * @param t Time
    * @return Exact x-velocity value
    */
  Real value_x(Real x, Real y, Real z, Real t) const
  {
    return std::sin(x * DX) * std::cos(y * DY) * std::sin(z * DZ) * std::sin(t);
  }

  /**
    * @brief Calculates exact solution for y-velocity component
    * @param x X-coordinate
    * @param y Y-coordinate
    * @param z Z-coordinate
    * @param t Time
    * @return Exact y-velocity value
    */
  Real value_y(Real x, Real y, Real z, Real t) const
  {
    return std::cos(x * DX) * std::sin(y * DY) * std::sin(z * DZ) * std::sin(t);
  }

  /**
    * @brief Calculates exact solution for z-velocity component
    * @param x X-coordinate
    * @param y Y-coordinate
    * @param z Z-coordinate
    * @param t Time
    * @return Exact z-velocity value
    */
  Real value_z(Real x, Real y, Real z, Real t) const
  {
    return 2 * std::cos(x * DX) * std::cos(y * DY) * std::cos(z * DZ) * std::sin(t);
  }

  /**
    * @brief Calculates exact solution for z-velocity component
    * @param x X-coordinate
    * @param y Y-coordinate
    * @param z Z-coordinate
    * @param t Time
    * @return Exact pressure value
    */
  Real value_p(Real x, Real y, Real z, Real t) const
  {
    return std::cos(x * DX) * std::cos(y * DY) * std::cos(z * DZ) * std::sin(t);
  }

};

/**
    * @brief Evaluates input string as mathematical expression
    * @param expr Mathematical expression as a string
    * @return Result of the mathematical expression
    */
inline double evaluateExpression(const std::string& expr) {
    std::string processedExpr = expr;
    size_t pos = processedExpr.find("M_PI");
    if (pos != std::string::npos) {
        std::ostringstream ss;
        ss << std::setprecision(16) << std::fixed << std::numbers::pi_v<long double>;
        processedExpr.replace(pos, 4, ss.str());
    }

    if (processedExpr.find('*') != std::string::npos) {
        std::istringstream iss(processedExpr);
        long double a, b;
        char op;
        iss >> a >> op >> b;
        return static_cast<double>(a * b);
    }

    return static_cast<double>(std::stold(processedExpr));
}


#endif
