#ifndef CONFIG_HPP
#define CONFIG_HPP

// Define the physical dimensions of the domain.
constexpr double LX = 1.0;
constexpr double LY = 1.0;
constexpr double LZ = 1.0;

// Define the number of cells in the x, y, z directions.
constexpr int NX = 100;
constexpr int NY = 100;
constexpr int NZ = 100;

// Define the space discretization in the x, y, z direction.
constexpr double DX = LX / NX;
constexpr double DY = LY / NY;
constexpr double DZ = LZ / NZ;

// Define the time step size.
constexpr double DT = 0.005;

// Define the Time interval.
constexpr double T = 1.0;

// Define the Reynolds number.
constexpr double RE = 400.0;

#endif // CONFIG_HPP
