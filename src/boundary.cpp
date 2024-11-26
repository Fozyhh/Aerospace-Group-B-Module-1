#include "boundary.hpp"

//<<<<<<< v2/array1D
// The method, which takes as input only the time step, is called by the program at the start of the time step.
// For each boundary node, it takes the exact value dor each component from the input and updates them.
// It also updates those values that are not directly on a face, but need an approximation.

