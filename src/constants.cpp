#include "constants.hpp"

Real SX = 0.0;
Real SY = 0.0;
Real SZ = 0.0;

// Domain lengths
Real LX = 10.0;  // Example value
Real LY = 20.0;  // Example value
Real LZ = 30.0;  // Example value

// Number of cells in each direction
int NX = 100, NY = 200, NZ = 300;

// Process grid dimensions
int PX = 4, PY = 4, PZ = 4;

// Space discretization steps
Real DX = LX / NX;
Real DY = LY / NY;
Real DZ = LZ / NZ;

// Time-related parameters
Real DT = 0.001;  // Example time step
Real Nt = 1000.0;   // Example total time
Real T = 0.1;  // Current time

// Boundary conditions
bool BX = false, BY = false, BZ = false; // Periodic by default

// Flow parameters
Real RE = 1000.0;  // Example Reynolds number
