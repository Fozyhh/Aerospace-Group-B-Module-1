# Aerospace-Group-B-Module-1
Repository for the project related to the "High Performance Scientific Computing in AeroSpace" course on Navier-Stokes equation for Incrompessible Fluids


# CPU Profiling, Memory Profiling and Memory check
In order to enable CPU profiling using GProf and Memory checks using Valgrind please follow these steps: 
Make sure before compiling that you have the `{PROJECT_ROOT/resources}` directory.

## Dependencies:
`sudo apt install binutils`
Valgrind:
`sudo apt install valgrind`
KCachegrind:
`sudo apt install kcachegrind`

## Normal build
``mkdir build
cd build
cmake ..
make
./main
``
## CPU profiling build
``mkdir build-profile
cd build-profile
cmake -DENABLE_PROFILING=ON ..
make profile
``
## Memory profiling build
`mkdir build-profiling
cd build-profiling
cmake -DENABLE_MEMPROF=ON ..
make profile`

## Memcheck build
``mkdir build-mem
cd build-mem
cmake -DENABLE_MEMCHECK=ON ..
make
ctest --verbose
``
You can then find the reports in the `resources` directory.

PS: It is not ideal to run both profiling and memcheck together as Gprof messes with the optimization flags and Valgrind greatly affects the execution speed of the executable.


# Error plot python script
## Dependencies
- matplotlib
- numpy

Run using: 
`` python3 vizualisation.py ``
